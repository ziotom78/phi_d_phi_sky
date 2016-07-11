{line 294 "phiskyandd.nw"}{ -*- mode: delphi -*- }
program CalcPhiD;

{$mode objfpc}{$h+}
uses Classes, SysUtils, INIFiles, DataTypes, Squeezer, Cfitsio, Healpix, Mpi,
     RotMatrix, SatelliteVelocities, ConvolvedParams;

{$linklib c} { Required by OpenMPI/MPICH }

const
    ProgramName = 'phid';
    
{line 1749 "phiskyandd.nw"}    MapColumns : Array[1..2] of Cfitsio.TColumn =
        ((Name: 'AVGPHI'; Count: 1; DataType: FitsTypeFloat; UnitStr: ''),
         (Name: 'HITS';   Count: 1; DataType: FitsTypeLong;  UnitStr: ''));

{line 307 "phiskyandd.nw"}type
    
{line 1403 "phiskyandd.nw"}TPercentage = 0..100;
TPercentageArray = array of TPercentage; { Open array }
{line 1410 "phiskyandd.nw"}TStringArray = array of String;  { Open array }
TDoubleArray = array of Double;  { Open array }
TBoolArray   = array of Boolean; { Open array }
{line 309 "phiskyandd.nw"}    
{line 461 "phiskyandd.nw"}TPhiDConfiguration = record
    PointingFileNames : TStringArray;
    SlMatrixFileNames : TStringArray;
    MbMatrixFileNames : TStringArray;
    SatelliteVelocityFileName : String;
    DipoleParams : TDipole;

    Quantiles : TPercentageArray;
    QuantileTableFileName : String;

    Nside : Uint16;
    OutputMapFileName : String;

    case SaveTods : Boolean of
    True: (TodFilePath : ShortString); { We need a ShortString here! }
end;
{line 719 "phiskyandd.nw"}TConvMatrixArray = array of TConvolutionMatrix;
TPhiDInputData = record
    MbMatrices : TConvMatrixArray;
    SlMatrices : TConvMatrixArray;
    SolSysVelEcl : TVector;
    SatelliteVelocity : TSatelliteVelocities;
end;
{line 941 "phiskyandd.nw"}TPhiDTod = record
    ObtTimes : Array of Int64;
    D : Array of Double;
    BslD : Array of Double;
    BmD : Array of Double;
    PhiD : Array of Double;
    Valid : Array of Boolean;
end;

{line 1797 "phiskyandd.nw"}procedure Log(const Message : String);
var
    MpiRank : Integer;

begin
    MpiRank := Mpi.CommRank(Mpi.World);
    WriteLn(
                Format('[%s #%d] %s',
                       [FormatDateTime('YYYY/MM/DD tt', Now),
                        MpiRank,
                        Message]));
end;
{line 1819 "phiskyandd.nw"}procedure EnsurePathExists(const FileName : String);
var
    Path : String;

begin
    Path := ExtractFilePath(ExpandFileName(FileName));
    if not ForceDirectories(Path) then
        raise Exception.CreateFmt('Unable to create the path "%s" ' +
                                  '(needed to save file "%s")',
                                  [Path, FileName]);
end;
{line 666 "phiskyandd.nw"}procedure PrintHelp;
begin
    if Mpi.CommRank(Mpi.World) = 0 then
        WriteLn(Format('Usage: %s PARAMETER_FILE', [ProgramName]));
end;
{line 858 "phiskyandd.nw"}procedure DivideMpiJobs(MpiRank, MpiSize : Integer;
                        NumOfFiles : Integer;
                        out FirstIdx, LastIdx : Integer);
begin
    FirstIdx := (NumOfFiles * MpiRank) div MpiSize;
    LastIdx := (NumOfFiles * (MpiRank + 1)) div MpiSize - 1;
end;
{line 1426 "phiskyandd.nw"}procedure InPlaceQuickSort(var A : Array of Double; FirstIdx, LastIdx : Integer);
Var
    i, j : LongInt;
    tmp, pivot : Double;

Begin
    i := FirstIdx;
    j := LastIdx;
    pivot := A[(FirstIdx + LastIdx) div 2];
    repeat
        while pivot > A[i] do Inc(i);
        While pivot < A[j] do Dec(j);
        if i <= j then begin
            tmp := A[i];
            A[i] := A[j];
            A[j] := tmp;
            Inc(i);
            Dec(j);
        end;
    until i > j;
    if FirstIdx < j then InPlaceQuickSort(A, FirstIdx, j);
    if i < LastIdx then InPlaceQuickSort(A, i, LastIdx);
End;
{line 1461 "phiskyandd.nw"}procedure AppendQuantiles(const A : TDoubleArray;
                          const Valid : TBoolArray;
                          const Quantiles : TPercentageArray;
                          var QuantileArray : TDoubleArray;
                          ChunkIdx : Integer);
var
    ValidValues : TDoubleArray;
    Idx, ValidIdx, QuantIdx : Integer;
    ValuesStr : String;
begin
    Assert(Length(A) = Length(Valid));
    Log(Format('Computing %d quantiles out of an array of %d elements…',
               [Length(Quantiles), Length(A)]));
    
{line 1485 "phiskyandd.nw"}SetLength(ValidValues, Length(A));
ValidIdx := Low(ValidValues);
for Idx := Low(A) to High(A) do
begin
    if Valid[Idx] then
    begin
        ValidValues[ValidIdx] := A[Idx];
        Inc(ValidIdx);
    end;
end;
if ValidIdx = Low(ValidValues) then
begin
    Log('…no valid values found for this OD, skipping the computation of quantiles');
    Exit;
end;

SetLength(ValidValues, ValidIdx);  { Truncate the tail of ValidValues }
Log(Format('…%d valid values found…', [Length(ValidValues)]));
{line 1475 "phiskyandd.nw"}    
{line 1513 "phiskyandd.nw"}                InPlaceQuickSort(ValidValues, Low(ValidValues), High(ValidValues));
Log(Format('…values have been sorted, their range is [%.3e, %.3e]…',
           [ValidValues[Low(ValidValues)], ValidValues[High(ValidValues)]]));
Idx := ChunkIdx * Length(Quantiles);
ValuesStr := '';
for QuantIdx := Low(Quantiles) to High(Quantiles) do
begin
    QuantileArray[Idx + QuantIdx] :=
        ValidValues[(Quantiles[QuantIdx] * Length(ValidValues)) div 100];
    if ValuesStr <> '' then ValuesStr := ValuesStr + ', ';
    ValuesStr := ValuesStr + Format('%.3e (%d%%)',
                                    [QuantileArray[Idx + QuantIdx], Quantiles[QuantIdx]]);
end;
{line 1476 "phiskyandd.nw"}    Log(Format('…the quantiles are %s.', [ValuesStr]));
end;
{line 1581 "phiskyandd.nw"}procedure ProjectTodOntoMap(const Pointings : TDetectorPointings;
                            const Tod : TDoubleArray;
                            const Valid : TBoolArray;
                            var BinnedMap : THealpixMap;
                            var HitMap : THealpixMap);
var
    Idx : Integer;
    PixelIdx : Cardinal;

begin
    Log(Format('Projecting %d samples into a NSIDE=%d map…',
               [Length(Tod), BinnedMap.Resolution.Nside]));
    for Idx := Low(Tod) to High(Tod) do
    begin
        if Valid[Idx] then
        begin
            PixelIdx := AnglesToPix(BinnedMap, Pointings.Theta[Idx], Pointings.Phi[Idx]);
            BinnedMap.Pixels[PixelIdx] := BinnedMap.Pixels[PixelIdx] + Tod[Idx];
            HitMap.Pixels[PixelIdx] := HitMap.Pixels[PixelIdx] + 1.0;
        end;
    end;
    Log('…projection completed.');
end;
{line 1630 "phiskyandd.nw"}procedure GatherQuantiles(const LocalQuantiles : Array of Double;
                          out OverallQuantiles : TDoubleArray);
var
    Idx : Integer;
    BufLengths : Array of Integer;
    Displacements : Array of Integer;
    LocalBufLength : Array[1..1] of Integer;

begin
    { Retrieve the number of quantiles computed by each MPI process }
    SetLength(BufLengths, Mpi.CommSize(Mpi.World));
    LocalBufLength[1] := Length(LocalQuantiles);
    Mpi.Gather(LocalBufLength, BufLengths, 0, Mpi.World);

    { Set up the array of displacements so that no holes will be left
      in OverallQuantiles }
    SetLength(Displacements, Length(BufLengths));
    Displacements[0] := 0;
    for Idx := Low(BufLengths) to High(BufLengths) - 1 do
        Displacements[Idx + 1] := Displacements[Idx] + BufLengths[Idx];

    { Gather the quantiles from each MPI process to the root process }
    SetLength(OverallQuantiles,
              Displacements[High(Displacements)] + BufLengths[High(BufLengths)]);
    Mpi.Gatherv(LocalQuantiles, OverallQuantiles, BufLengths,
                Displacements, 0, Mpi.World);
end;
{line 1660 "phiskyandd.nw"}procedure SaveQuantileTOD(const FileName : String;
                          const Quantiles : TPercentageArray;
                          const QuantileArray : TDoubleArray);
const
    PercentTableDef : Array[1..1] of Cfitsio.TColumn =
        ((Name: 'PERCENT'; Count: 1; DataType: FitsTypeShort; UnitStr: 'Percentage'));

var
    QuantileTableDef : Array of Cfitsio.TColumn;
    QuantIdx, Idx : Integer;
    F : TFitsFile;
    CurQuantiles : TDoubleArray;

begin
    SetLength(QuantileTableDef, Length(Quantiles));
    for QuantIdx := Low(QuantileTableDef) to High(QuantileTableDef) do
    begin
        with QuantileTableDef[QuantIdx] do
        begin
            Name := Format('Q%.3d', [Quantiles[QuantIdx]]);
            Count := 1;
            DataType := FitsTypeFloat;
            UnitStr := '';
        end;
    end;

    try
        Log(Format('Saving quantiles into FITS file "%s"…', [FileName]));
        F := Cfitsio.CreateFile(FileName, Overwrite);
        try
            CreateTable(F, BinaryTable, 0, PercentTableDef, 'PERCENTAGES');
            WriteColumn(F, 1, 1, 1, Quantiles);

            CreateTable(F, BinaryTable, 0, QuantileTableDef, 'QUANTILES');
            SetLength(CurQuantiles, Length(QuantileArray) div Length(Quantiles));
            for QuantIdx := 0 to Length(Quantiles) - 1 do
            begin
                for Idx := 0 to Length(CurQuantiles) - 1 do
                    CurQuantiles[Idx] := QuantileArray[Idx * Length(Quantiles) + QuantIdx];
                WriteColumn(F, 1 + QuantIdx, 1, 1, CurQuantiles);
            end;
        finally
            Cfitsio.CloseFile(F);
        end;
    except
        on E : Exception do Log(Format('Unable to write file "%s": %s',
                                       [FileName, E.Message]));
    end;
    Log(Format('…file "%s" saved.', [FileName]));
end;
{line 1849 "phiskyandd.nw"}procedure GetAllValuesFromSection(IniFile : TIniFile;
                                  SectionName : String;
                                  var ValueList : TStringArray);
var
    KeyList  : TStringList;
    Idx      : Integer;

begin
    KeyList := TStringList.Create;
    try
        IniFile.ReadSection(SectionName, KeyList);
        SetLength(ValueList, KeyList.Count);
        for Idx := 0 to KeyList.Count - 1 do
            ValueList[Idx] := (IniFile.ReadString(SectionName,
                                                  KeyList[Idx], ''));
    finally
        KeyList.Free;
    end;
end;
{line 1876 "phiskyandd.nw"}procedure GetSequenceOfFiles(IniFile : TIniFile;
                             const SectName : String;
                             var FileNames : TStringArray);
var
    FirstIdx, LastIdx, Idx : Integer;
    Template : String;

begin
    if IniFile.ValueExists(SectName, 'first_index') and
       IniFile.ValueExists(SectName, 'last_index') and
       IniFile.ValueExists(SectName, 'template') then
    begin
        FirstIdx := IniFile.ReadInteger(SectName, 'first_index', 0);
        LastIdx := IniFile.ReadInteger(SectName, 'last_index', 0);
        Template := IniFile.ReadString(SectName, 'template', '');

        SetLength(FileNames, LastIdx - FirstIdx + 1);
        for Idx := FirstIdx to LastIdx do
            FileNames[Idx - FirstIdx] := Format(Template, [Idx]);
    end else
        GetAllValuesFromSection(IniFile, SectName, FileNames);
end;
{line 1905 "phiskyandd.nw"}procedure GetPercentages(const InputStr : String;
                         out Percentages : TPercentageArray);
var
    PercStrList : TStringList;
    Code : Word;
    Idx : Integer;

begin
    PercStrList := TStringList.Create;
    try
        ExtractStrings([','], [' ', #9], PChar(InputStr), PercStrList);
        SetLength(Percentages, PercStrList.Count);
        for Idx := 0 to PercStrList.Count - 1 do
        begin
            Val(PercStrList[Idx], Percentages[Idx], Code);
            if Code <> 0 then
                raise ERangeError(Format('"%s" is not a percentage',
                                         [PercStrList[Idx]]));
        end;
    finally
        PercStrList.Free;
    end;
end;
{line 2039 "phiskyandd.nw"}function ReadBoolean(const IniFile : TIniFile;
                     const SectName : String;
                     const KeyName : String;
                     DefaultValue : Boolean) : Boolean;
var
    DefaultValStr : String;
    Value : String;

begin
    if DefaultValue then DefaultValStr := 'true' else DefaultValStr := 'false';

    Value := UpperCase(IniFile.ReadString(SectName, KeyName,
                                          DefaultValStr));
    if (Value = 'TRUE') or (Value = 'YES') or (Value = 'ON') then
        Exit(True)
    else if (Value = 'FALSE') or (Value = 'NO') or (Value = 'OFF') then
        Exit(False);

    { Unable to understand what the user wants, so stick with the default }
    Exit(DefaultValue);
end;
{line 2121 "phiskyandd.nw"}procedure WriteFileNames(const SectName : String;
                         const NameList : TStringArray);
var
    Idx : Integer;

begin
    WriteLn(SectName);
    for Idx := 0 to Length(NameList) - 1 do
        WriteLn(Format('file_name_%.4d = %s',
                       [Idx, NameList[Idx]]));

    WriteLn;
end;
{line 2142 "phiskyandd.nw"}procedure WriteQuantiles(const KeyName : String;
                         const Quantiles : TPercentageArray);
var
    Idx : Integer;

begin
    Write(KeyName, ' = ');
    for Idx := 0 to Length(Quantiles) - 1 do
    begin
        Write(Quantiles[Idx]);
        if Idx < High(Quantiles) then
            Write(', ')
        else
            WriteLn;
    end;
end;
{line 982 "phiskyandd.nw"}procedure ComputeD(const SolSysVelEcl : TVector;
                   const Pointings : TDetectorPointings;
                   const Vel : TSatelliteVelocities;
                   var D : TDoubleArray);
var
    Idx : Integer;

begin
    SetLength(D, Length(Pointings.Theta));
    for Idx := 0 to Length(Pointings.Theta) - 1 do
    begin
        with Pointings do
            D[Idx] := Convolve(DiracDelta, SolSysVelEcl, Vel,
                               ScetTimes[Idx], Theta[Idx], Phi[Idx], Psi[Idx]);
    end;
end;
{line 1005 "phiskyandd.nw"}procedure CalculatePhiD(var PhiDTod : TPhiDTod);
var
    Idx : Integer;
    BmDDiff : Double;

begin
    with PhiDTod do
    begin
        SetLength(PhiD, Length(BslD));
        SetLength(Valid, Length(BmD));
        for Idx := 1 to Length(PhiD) - 1 do
        begin
            BmDDiff := BmD[Idx] - BmD[Idx - 1];
            if BmDDiff <> 0.0 then
            begin
                PhiD[Idx] := (BslD[Idx] - BslD[Idx - 1]) / BmDDiff;
                Valid[Idx] := True;
            end else
                Valid[Idx] := False;
        end;
    end;
end;
{line 1068 "phiskyandd.nw"}{ Forward declaration }
procedure SumConvMatricesIntensities(const Matrices : TConvMatrixArray;
                                     const SolSysVelEcl : TVector;
                                     const SatVel : TSatelliteVelocities;
                                     const Pointings : TDetectorPointings;
                                     var DestVector : TDoubleArray); forward;

procedure ProcessPointingFile(const FileName : String;
                              const InputData : TPhiDInputData;
                              out PhiDTod : TPhiDTod;
                              var BinnedMap : THealpixMap;
                              var HitMap : THealpixMap);
var
    FileHeader : Squeezer.TFileHeader;
    Pointings : TDetectorPointings;

begin
    Log(Format('Reading pointing file %s…', [FileName]));
    ReadDetectorPointings(FileName, FileHeader, Pointings);
    Log(Format('…pointing file read, %d samples found',
               [Length(Pointings.ObtTimes)]));

    SetLength(PhiDTod.ObtTimes, Length(Pointings.ObtTimes));
    Move(Pointings.ObtTimes[0], PhiDTod.ObtTimes[0],
         Length(Pointings.ObtTimes) * SizeOf(Pointings.ObtTimes[0]));

    Log('Computing the term D…');
    ComputeD(InputData.SolSysVelEcl, Pointings,
             InputData.SatelliteVelocity, PhiDTod.D);
    Log('…term D computed');

    Log(Format('Calculating convolutions using %d+%d (MB/SL) convolution matrices…',
              [Length(InputData.MbMatrices), Length(InputData.SlMatrices)]));
    SumConvMatricesIntensities(InputData.MbMatrices, InputData.SolSysVelEcl,
                               InputData.SatelliteVelocity, Pointings, PhiDTod.BmD);
    SumConvMatricesIntensities(InputData.SlMatrices, InputData.SolSysVelEcl,
                               InputData.SatelliteVelocity, Pointings, PhiDTod.BslD);
    Log('…convolutions calculated');

    Log('Calculating φ_D…');
    CalculatePhiD(PhiDTod);
    Log('…φ_D calculated');

    ProjectTodOntoMap(Pointings, PhiDTod.PhiD, PhiDTod.Valid, BinnedMap, HitMap);
end;
{line 1293 "phiskyandd.nw"}procedure SavePhiTod(const FileName : String;
                     const PhiDTod : TPhiDTod);
const
    FileColumns : Array[1..5] of Cfitsio.TColumn =
        ((Name: 'OBTTIME'; Count: 1; DataType: FitsTypeDouble;  UnitStr: ''),
         (Name: 'BSLD';    Count: 1; DataType: FitsTypeFloat;   UnitStr: 'K_CMB'),
         (Name: 'BMD';     Count: 1; DataType: FitsTypeFloat;   UnitStr: 'K_CMB'),
         (Name: 'PHID';    Count: 1; DataType: FitsTypeFloat;   UnitStr: 'K_CMB'),
         (Name: 'VALID';   Count: 1; DataType: FitsTypeLogical; UnitStr: ''));

var
    F : TFitsFile;

begin
    Log(Format('Saving %d values of φ_D into file %s…',
               [Length(PhiDTod.ObtTimes), FileName]));
    try
        EnsurePathExists(FileName);
        F := Cfitsio.CreateFile(FileName, Overwrite);
        try
            Cfitsio.CreateTable(F, BinaryTable, 0, FileColumns, 'PHID');
            Cfitsio.WriteComment(F, 'File created by the phid program');
            Cfitsio.WriteComment(F, 'Phid was compiled on ' +
                                        {$i %date} + ' ' + {$i %time});

            Cfitsio.WriteColumn(F, 1, 1, 1, PhiDTod.ObtTimes);
            Cfitsio.WriteColumn(F, 2, 1, 1, PhiDTod.BslD);
            Cfitsio.WriteColumn(F, 3, 1, 1, PhiDTod.BmD);
            Cfitsio.WriteColumn(F, 4, 1, 1, PhiDTod.PhiD);
            Cfitsio.WriteColumn(F, 5, 1, 1, PhiDTod.Valid);

            Log(Format('…done, file %s has been saved', [FileName]));
        finally
            Cfitsio.CloseFile(F);
        end;
    except
        on E : EFitsError do Log(Format('Unable to write file %s: %s',
                                        [FileName, E.message]));
    end;
end;
{line 1935 "phiskyandd.nw"}procedure ReadConfiguration(const FileName : String;
                            out Configuration : TPhiDConfiguration);
var
    IniFile : TIniFile;
    ListOfPercentages : String;
begin
    IniFile := TIniFile.Create(FileName);
    try
        with Configuration do
        begin
            GetAllValuesFromSection(IniFile, 'Main beam matrices', MbMatrixFileNames);
            GetAllValuesFromSection(IniFile, 'Sidelobe matrices', SlMatrixFileNames);

            GetSequenceOfFiles(IniFile, 'Pointings', PointingFileNames);

            SatelliteVelocityFileName :=
                IniFile.ReadString('Input', 'satellite_velocity_file', '');
            DipoleParams.DirTheta :=
                IniFile.ReadFloat('Input', 'dipole_dir_theta_ecl', 1.7656131194951572);
            DipoleParams.DirPhi :=
                IniFile.ReadFloat('Input', 'dipole_dir_phi_ecl', 2.9958896005735780);
            DipoleParams.SpeedMS :=
                IniFile.ReadFloat('Input', 'dipole_speed_m_s', 370082.2332);

            ListOfPercentages :=
                IniFile.ReadString('Output', 'quantiles', '25,50,75');
            GetPercentages(ListOfPercentages, Quantiles);
            QuantileTableFileName := IniFile.ReadString('Output', 'quantiles_file_name', '');
            Nside := IniFile.ReadInteger('Output', 'map_nside', 64);
            if not Healpix.IsNsideValid(Nside) then
                raise Exception.CreateFmt('Invalid NSIDE value (%d)',
                                          [Nside]);

            OutputMapFileName :=
                IniFile.ReadString('Output', 'output_file_name', 'phi_d.fits');

            SaveTods := ReadBoolean(IniFile, 'Output', 'save_tods', False);
            if SaveTods then
                TodFilePath :=
                    IniFile.ReadString('Output', 'tod_file_path', './');
        end;
    finally
        IniFile.Free;
    end;
end;
{line 2068 "phiskyandd.nw"}procedure PrintConfiguration(const Configuration : TPhiDConfiguration);
begin
    with Configuration do
    begin
        WriteFileNames('[Pointings]', PointingFileNames);
        WriteFileNames('[Main beam matrices]', MbMatrixFileNames);
        WriteFileNames('[Sidelobe matrices]', SlMatrixFileNames);

        WriteLn('[Input]');
        WriteLn('satellite_velocity = ', SatelliteVelocityFileName);
        WriteLn;

        WriteLn('[Output]');
        WriteLn('nside = ', Nside);
        WriteLn('output_file_name = ', OutputMapFileName);
        WriteQuantiles('quantiles', Quantiles);
        WriteLn('quantiles_file_name', QuantileTableFileName);
        WriteLn('save_tods = ', SaveTods);
        if SaveTods then
            WriteLn('tod_file_path = ', TodFilePath);
    end;
end;
{line 2168 "phiskyandd.nw"}procedure LoadConvolutionMatrices(const FileNames : TStringArray;
                                  out Matrices : TConvMatrixArray);
var
    Idx : Integer;
    CurFileName : String;

begin
    SetLength(Matrices, Length(FileNames));
    for Idx := 0 to Length(FileNames) - 1 do
    begin
        CurFileName := FileNames[Idx];
        Log(Format('Reading file %s', [CurFileName]));
        LoadConvolutionMatrix(CurFileName, Matrices[Idx]);
    end;
end;
{line 2186 "phiskyandd.nw"}procedure SumConvMatricesIntensities(const Matrices : TConvMatrixArray;
                                     const SolSysVelEcl : TVector;
                                     const SatVel : TSatelliteVelocities;
                                     const Pointings : TDetectorPointings;
                                     var DestVector : TDoubleArray);
var
    MatIdx, Idx : Integer;

begin
    SetLength(DestVector, Length(Pointings.Theta));
    for Idx := 0 to Length(DestVector) - 1 do
        DestVector[Idx] := 0.0;

    for MatIdx := 0 to Length(Matrices) - 1 do
    begin
        Log(Format('Applying convolution matrix %d/%d', [Idx + 1, Length(Matrices)]));
        for Idx := 0 to Length(Pointings.Theta) - 1 do
            DestVector[Idx] := DestVector[Idx] +
                Convolve(Matrices[MatIdx], SolSysVelEcl, SatVel,
                         Pointings.ObtTimes[Idx],
                         Pointings.Theta[Idx], Pointings.Phi[Idx], Pointings.Psi[Idx]);
    end;
end;

{line 315 "phiskyandd.nw"}var
    
{line 878 "phiskyandd.nw"}FirstFileIdx, LastFileIdx : Integer;
MpiRank, MpiSize : Integer;
Idx : Integer;
{line 1206 "phiskyandd.nw"}QuantileArray, OverallQuantileArray : TDoubleArray;
{line 1281 "phiskyandd.nw"}CurFileIdx : Integer;
{line 1560 "phiskyandd.nw"}BinnedMap, HitMap : THealpixMap;
OverallBinnedMap, OverallHitMap : THealpixMap;
{line 1786 "phiskyandd.nw"}FitsFile : TFitsFile;
{line 317 "phiskyandd.nw"}    
{line 689 "phiskyandd.nw"}Configuration : TPhiDConfiguration;
{line 739 "phiskyandd.nw"}InputData : TPhiDInputData;
{line 971 "phiskyandd.nw"}PhiTod : TPhiDTod;

{line 319 "phiskyandd.nw"}begin
    Mpi.Init;
    try
        try
            
{line 651 "phiskyandd.nw"}if ParamCount <> 1 then
begin
   PrintHelp;
   Exit;
end;
{line 678 "phiskyandd.nw"}ReadConfiguration(ParamStr(1), Configuration);
if Mpi.CommRank(Mpi.World) = 0 then
    PrintConfiguration(Configuration);
{line 891 "phiskyandd.nw"}MpiRank := Mpi.CommRank(Mpi.World);
MpiSize := Mpi.CommSize(Mpi.World);
if MpiSize > Length(Configuration.PointingFileNames) then
begin
    MpiSize := Length(Configuration.PointingFileNames);
    if MpiRank >= MpiSize then
    begin
        Log(Format('Too many MPI processes (%d) and too few ' +
                   'files (%d): this process will stop',
                   [MpiSize, Length(Configuration.PointingFileNames)]));
        Exit;
    end;
end;

DivideMpiJobs(MpiRank, MpiSize, Length(Configuration.PointingFileNames),
              FirstFileIdx, LastFileIdx);

case LastFileIdx - FirstFileIdx + 1 of
0: Log('No files to load');
1: Log(Format('This MPI process will analyze pointing file %s',
              [Configuration.PointingFileNames[FirstFileIdx]]));
else Log(Format('This MPI process will analyze %d pointing files ' +
                    '(out of %d): %s … %s',
                [LastFileIdx - FirstFileIdx + 1,
                 Length(Configuration.PointingFileNames),
                 Configuration.PointingFileNames[FirstFileIdx],
                 Configuration.PointingFileNames[LastFileIdx]]))
end;
{line 780 "phiskyandd.nw"}with Configuration do
begin
    Log('Loading the convolution matrices');

    LoadConvolutionMatrices(MbMatrixFileNames, InputData.MbMatrices);
    LoadConvolutionMatrices(SlMatrixFileNames, InputData.SlMatrices);

    Log('Convolution matrices loaded, loading the satellite velocities');

    LoadFromFile(SatelliteVelocityFileName, InputData.SatelliteVelocity);

    Log('Satellite velocities loaded');

    { Set up the velocity of the satellite wrt the Solar System }
    with InputData.SolSysVelEcl do
        Healpix.AnglesToVector(DipoleParams.DirTheta, DipoleParams.DirPhi, x, y, z);
    InputData.SolSysVelEcl :=
        ScaleVector(InputData.SolSysVelEcl, DipoleParams.SpeedMS);

end;
{line 1213 "phiskyandd.nw"}SetLength(QuantileArray, (LastFileIdx - FirstFileIdx + 1) * Length(Configuration.Quantiles));
{line 1570 "phiskyandd.nw"}InitMap(Configuration.Nside, Ring, BinnedMap);
InitMap(Configuration.Nside, Ring, HitMap);

InitMap(Configuration.Nside, Ring, OverallBinnedMap);
InitMap(Configuration.Nside, Ring, OverallHitMap);
{line 1219 "phiskyandd.nw"}for CurFileIdx := FirstFileIdx to LastFileIdx do
begin
    with Configuration do
    begin
        try
            ProcessPointingFile(PointingFileNames[CurFileIdx],
                                InputData, PhiTod, BinnedMap, HitMap);

            AppendQuantiles(PhiTod.PhiD, PhiTod.Valid, Configuration.Quantiles,
                            QuantileArray, CurFileIdx - FirstFileIdx);

            if Configuration.SaveTods then
                SavePhiTod(ConcatPaths([TodFilePath,
                    ChangeFileExt(ExtractFileName(PointingFileNames[CurFileIdx]),
                                  '-phid.fits')]),
                              PhiTod);
        except
            on E : Exception do
                Log(Format('Unable to process file "%s" (%s), skipping…',
                           [PointingFileNames[CurFileIdx],
                            E.Message]));
        end;
    end;
end;
{line 1713 "phiskyandd.nw"}GatherQuantiles(QuantileArray, OverallQuantileArray);
if MpiRank = 0 then
begin
    EnsurePathExists(Configuration.QuantileTableFileName);
    SaveQuantileTOD(Configuration.QuantileTableFileName,
                    Configuration.Quantiles, OverallQuantileArray);
end;
{line 1729 "phiskyandd.nw"}Log(Format('Reducing %d binned values…', [Length(BinnedMap.Pixels)]));
Mpi.Reduce(BinnedMap.Pixels, OverallBinnedMap.Pixels, MPI_SUM, 0, Mpi.World);
Log(Format('Reducing %d hit count pixels…', [Length(BinnedMap.Pixels)]));
Mpi.Reduce(HitMap.Pixels, OverallHitMap.Pixels, MPI_SUM, 0, Mpi.World);

if MpiRank = 0 then
begin
    for Idx := 0 to Length(BinnedMap.Pixels) - 1 do
    begin
        if OverallHitMap.Pixels[Idx] > 0 then
            OverallBinnedMap.Pixels[Idx] :=
                OverallBinnedMap.Pixels[Idx] / OverallHitMap.Pixels[Idx]
        else
            OverallBinnedMap.Pixels[Idx] := Healpix.Unseen;
    end;
end;
{line 1762 "phiskyandd.nw"}if MpiRank = 0 then
begin
    Log(Format('Writing binned map in %s…', [Configuration.OutputMapFileName]));
    try
        EnsurePathExists(Configuration.OutputMapFileName);
        FitsFile := Cfitsio.CreateFile(Configuration.OutputMapFileName, Overwrite);
        try
            Cfitsio.CreateTable(FitsFile, BinaryTable, 0, MapColumns, 'PHIMAP');
            Healpix.WriteHealpixKeywords(FitsFile, OverallBinnedMap);
            Cfitsio.WriteColumn(FitsFile, 1, 1, 1, OverallBinnedMap.Pixels);
            Cfitsio.WriteColumn(FitsFile, 2, 1, 1, OverallHitMap.Pixels);
        finally
            Cfitsio.CloseFile(FitsFile);
        end;
        Log('…done, map has been written successfully');
    except
        on E : Exception do Log('…error, unable to write the map: ' + E.Message);
    end;
end;
{line 324 "phiskyandd.nw"}        finally
            Mpi.Finalize;
        end;
    except
        on E : Exception do WriteLn('Error: ' + E.message);
    end;
end.
