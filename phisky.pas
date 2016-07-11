{line 339 "phiskyandd.nw"}{ -*- mode: delphi -*- }
program CalcPhiSky;

{$mode objfpc}{$h+}
uses Classes, SysUtils, INIFiles, DataTypes, Ringsets, Squeezer, Cfitsio, Healpix, Mpi;

{$linklib c} { Required by OpenMPI/MPICH }

const
    ProgramName = 'phisky';
    
{line 1749 "phiskyandd.nw"}    MapColumns : Array[1..2] of Cfitsio.TColumn =
        ((Name: 'AVGPHI'; Count: 1; DataType: FitsTypeFloat; UnitStr: ''),
         (Name: 'HITS';   Count: 1; DataType: FitsTypeLong;  UnitStr: ''));

{line 351 "phiskyandd.nw"}type
    
{line 1403 "phiskyandd.nw"}TPercentage = 0..100;
TPercentageArray = array of TPercentage; { Open array }
{line 1410 "phiskyandd.nw"}TStringArray = array of String;  { Open array }
TDoubleArray = array of Double;  { Open array }
TBoolArray   = array of Boolean; { Open array }
{line 353 "phiskyandd.nw"}    
{line 483 "phiskyandd.nw"}TPhiSkyConfiguration = record
    MbRingsetFileNames : TStringArray;
    SlRingsetFileNames : TStringArray;
    PointingFileNames : TStringArray;
    TemperatureFileNames : TStringArray;

    QualityFlagMask : UInt32;

    InterpolationOrder : Byte;
    Quantiles : TPercentageArray;
    QuantileTableFileName : String;

    Nside : Uint16;
    OutputMapFileName : String;

    case SaveTods : Boolean of
    True: (TodFilePath : ShortString); { We need a ShortString here! }
end;
{line 729 "phiskyandd.nw"}TRingsetArray = array of TRingset;
TPhiSkyInputData = record
    MbRingsets : TRingsetArray;
    SlRingsets : TRingsetArray;
end;
{line 958 "phiskyandd.nw"}TPhiSkyTod = record
    ObtTimes : Array of Int64;
    TskyMeas : Array of Double;
    BslTsky : Array of Double;
    BmTsky : Array of Double;
    PhiSky : Array of Double;
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
{line 1032 "phiskyandd.nw"}procedure CalculatePhiSky(const Flags : Array of UInt32;
                          QualityFlagMask : UInt32;
                          var PhiSkyTod : TPhiSkyTod);
var
    Idx : Integer;

begin
    with PhiSkyTod do
    begin
        SetLength(PhiSky, Length(BslTsky));
        SetLength(Valid, Length(BslTsky));
        for Idx := 0 to Length(PhiSky) - 1 do
        begin
            if (TskyMeas[Idx] <> 0.0) and ((Flags[Idx] and QualityFlagMask) = 0) then
            begin
                PhiSky[Idx] := BslTsky[Idx] / TskyMeas[Idx];
                Valid[Idx] := True;
            end else
                Valid[Idx] := False;
        end;
    end;
end;
{line 1124 "phiskyandd.nw"}{ Forward declaration }
procedure SumRingsetIntensities(const Ringsets : TRingsetArray;
                                const Pointings : TDetectorPointings;
                                var DestVector : TDoubleArray); forward;

procedure ProcessFilePair(const PntFileName, TskyFileName : String;
                          const InputData : TPhiSkyInputData;
                          QualityFlagMask : UInt32;
                          out PhiSkyTod : TPhiSkyTod;
                          var BinnedMap : THealpixMap;
                          var HitMap : THealpixMap);

var
    FileHeader : Squeezer.TFileHeader;
    Pointings : TDetectorPointings;
    Temperature : TDifferencedData;

begin
    Log(Format('Reading pointing file %s…', [PntFileName]));
    ReadDetectorPointings(PntFileName, FileHeader, Pointings);
    Log(Format('…pointing file read, %d samples found',
               [Length(Pointings.ObtTimes)]));

    Log(Format('Reading temperature file %s…', [TskyFileName]));
    ReadDifferencedData(TskyFileName, FileHeader, Temperature);
    Log(Format('…temperature file read, %d samples found',
               [Length(Pointings.ObtTimes)]));

    SetLength(PhiSkyTod.ObtTimes, Length(Pointings.ObtTimes));
    Move(Pointings.ObtTimes[0], PhiSkyTod.ObtTimes[0],
         Length(Pointings.ObtTimes) * SizeOf(Pointings.ObtTimes[0]));

    SetLength(PhiSkyTod.TskyMeas, Length(Temperature.SkyLoad));
    Move(Temperature.SkyLoad[0], PhiSkyTod.TskyMeas[0],
         Length(Temperature.SkyLoad) * SizeOf(Temperature.SkyLoad[0]));

    Log(Format('Calculating convolutions using %d+%d (MB/SL) ringsets…',
              [Length(InputData.MbRingsets), Length(InputData.SlRingsets)]));
    SumRingsetIntensities(InputData.MbRingsets, Pointings, PhiSkyTod.BmTsky);
    SumRingsetIntensities(InputData.SlRingsets, Pointings, PhiSkyTod.BslTsky);
    Log('…convolutions calculated');

    Log('Calculating φ_sky…');
    CalculatePhiSky(Temperature.Flags, QualityFlagMask, PhiSkyTod);
    Log('…φ_sky calculated');

    ProjectTodOntoMap(Pointings, PhiSkyTod.PhiSky, PhiSkyTod.Valid, BinnedMap, HitMap);
end;
{line 1341 "phiskyandd.nw"}procedure SavePhiTod(const FileName : String;
                     const PhiSkyTod : TPhiSkyTod);
const
    FileColumns : Array[1..6] of Cfitsio.TColumn =
        ((Name: 'OBTTIME'; Count: 1; DataType: FitsTypeDouble;  UnitStr: ''),
         (Name: 'SLCONV';  Count: 1; DataType: FitsTypeFloat;   UnitStr: 'K_CMB'),
         (Name: 'MBCONV';  Count: 1; DataType: FitsTypeFloat;   UnitStr: 'K_CMB'),
         (Name: 'TSKY';    Count: 1; DataType: FitsTypeFloat;   UnitStr: 'K_CMB'),
         (Name: 'PHISKY';  Count: 1; DataType: FitsTypeFloat;   UnitStr: 'K_CMB'),
         (Name: 'VALID';   Count: 1; DataType: FitsTypeLogical; UnitStr: ''));

var
    F : TFitsFile;

begin
    Log(Format('Saving %d values of φ_sky into file %s',
               [Length(PhiSkyTod.ObtTimes), FileName]));
    try
        EnsurePathExists(FileName);
        F := Cfitsio.CreateFile(FileName, Overwrite);
        try
            Cfitsio.CreateTable(F, BinaryTable, 0, FileColumns, 'PHISKY');
            Cfitsio.WriteComment(F, 'File created by the phisky program');
            Cfitsio.WriteComment(F, 'Phisky was compiled on ' +
                                        {$i %date} + ' ' + {$i %time});
            Cfitsio.WriteColumn(F, 1, 1, 1, PhiSkyTod.ObtTimes);
            Cfitsio.WriteColumn(F, 2, 1, 1, PhiSkyTod.BslTsky);
            Cfitsio.WriteColumn(F, 3, 1, 1, PhiSkyTod.BmTsky);
            Cfitsio.WriteColumn(F, 4, 1, 1, PhiSkyTod.TskyMeas);
            Cfitsio.WriteColumn(F, 5, 1, 1, PhiSkyTod.PhiSky);
            Cfitsio.WriteColumn(F, 6, 1, 1, PhiSkyTod.Valid);
        finally
            Cfitsio.CloseFile(F);
        end;
    except
        on E : EFitsError do Log(Format('Unable to write file %s: %s',
                                        [FileName, E.message]));
    end;
end;
{line 1985 "phiskyandd.nw"}procedure ReadConfiguration(const FileName : String;
                            out Configuration : TPhiSkyConfiguration);
var
    IniFile : TIniFile;
    ListOfPercentages : String;
begin
    IniFile := TIniFile.Create(FileName);
    try
        with Configuration do
        begin
            GetAllValuesFromSection(IniFile, 'Main ringsets', MbRingsetFileNames);
            GetAllValuesFromSection(IniFile, 'Side ringsets', SlRingsetFileNames);

            GetSequenceOfFiles(IniFile, 'Pointings', PointingFileNames);
            GetSequenceOfFiles(IniFile, 'Temperatures', TemperatureFileNames);
            if Length(PointingFileNames) <> Length(TemperatureFileNames) then
                raise Exception.CreateFmt('# of pointing/temperature files differ (%d/%d)',
                                          [Length(PointingFileNames),
                                           Length(TemperatureFileNames)]);

            QualityFlagMask := IniFile.ReadInteger('Input', 'quality_flag_mask', 6111248);

            InterpolationOrder :=
                IniFile.ReadInteger('Output', 'interpolation_order', 5);
            ListOfPercentages :=
                IniFile.ReadString('Output', 'quantiles', '25,50,75');
            GetPercentages(ListOfPercentages, Quantiles);
            QuantileTableFileName := IniFile.ReadString('Output', 'quantiles_file_name', '');
            Nside := IniFile.ReadInteger('Output', 'map_nside', 64);
            if not Healpix.IsNsideValid(Nside) then
                raise Exception.CreateFmt('Invalid NSIDE value (%d)',
                                          [Nside]);

            OutputMapFileName :=
                IniFile.ReadString('Output', 'output_file_name', 'phi_sky.fits');

            SaveTods := ReadBoolean(IniFile, 'Output', 'save_tods', False);
            if SaveTods then
                TodFilePath :=
                    IniFile.ReadString('Output', 'tod_file_path', './');
        end;
    finally
        IniFile.Free;
    end;
end;
{line 2093 "phiskyandd.nw"}procedure PrintConfiguration(const Configuration : TPhiSkyConfiguration);
begin
    with Configuration do
    begin
        WriteFileNames('[Pointings]', PointingFileNames);
        WriteFileNames('[Temperatures]', TemperatureFileNames);
        WriteFileNames('[Main ringsets]', MbRingsetFileNames);
        WriteFileNames('[Side ringsets]', SlRingsetFileNames);

        WriteLn('[Output]');
        WriteLn('interpolation_order = ', InterpolationOrder);
        WriteLn('nside = ', Nside);
        WriteLn('output_file_name = ', OutputMapFileName);
        WriteQuantiles('quantiles', Quantiles);
        WriteLn('quantiles_file_name', QuantileTableFileName);
        WriteLn('save_tods = ', SaveTods);
        if SaveTods then
            WriteLn('tod_file_path = ', TodFilePath);
    end;
end;
{line 2218 "phiskyandd.nw"}procedure LoadRingsets(const FileNames : TStringArray;
                       InterpolationOrder : Integer;
                       out Ringsets : TRingsetArray);
var
    Idx : Integer;
    CurFileName : String;

begin
    SetLength(Ringsets, Length(FileNames));
    for Idx := 0 to Length(FileNames) - 1 do
    begin
        CurFileName := FileNames[Idx];
        Log(Format('Reading file %s', [CurFileName]));
        LoadRingsetFromFile(CurFileName,
                            InterpolationOrder,
                            Ringsets[Idx]);
    end;
end;
{line 2246 "phiskyandd.nw"}procedure UpdateUserAboutStatus(Idx, NumOfElements : Integer; Data : Pointer);
var
    Percentage : Real;

begin
    if NumOfElements > 0 then
        Percentage := (Idx * 100.0) / (NumOfElements - 1)
    else
        Percentage := 100.0;

    Log(Format('   Ringset application: %d/%d (%.1f%%)',
               [Idx + 1, NumOfElements, Percentage]));
end;

procedure SumRingsetIntensities(const Ringsets : TRingsetArray;
                                const Pointings : TDetectorPointings;
                                var DestVector : TDoubleArray);
const
    LogUpdateDelayInMs = 30000.0; { = 30 s }
var
    Idx : Integer;
begin
    SetLength(DestVector, Length(Pointings.Theta));
    for Idx := 0 to Length(DestVector) - 1 do
        DestVector[Idx] := 0.0;

    for Idx := 0 to Length(Ringsets) - 1 do
    begin
        Log(Format('Applying ringset %d/%d', [Idx + 1, Length(Ringsets)]));
        GetRingsetIntensities(Ringsets[Idx], Pointings.Theta,
                              Pointings.Phi, Pointings.Psi,
                              DestVector, AddElements,
                              @UpdateUserAboutStatus, LogUpdateDelayInMs, nil);
    end;
end;
{line 2292 "phiskyandd.nw"}function MemoryUsed(const InputData : TPhiSkyInputData) : String;
const
    Units : Array[1..4] of ShortString = ('bytes', 'KB', 'MB', 'GB');

var
    Size : Real = 0.0;
    FormattedNum : String;
    Idx, UnitIdx : Integer;

begin
    with InputData do
    begin
        for Idx := 0 to Length(MbRingsets) - 1 do
            Size := Size + BytesUsed(MbRingsets[Idx]);

        for Idx := 0 to Length(SlRingsets) -  1 do
            Size := Size + BytesUsed(SlRingsets[Idx]);
    end;

    UnitIdx := Low(Units);
    while (Size > 10 * 1024) and (UnitIdx < High(Units)) do
    begin
        Size := Size / 1024;
        Inc(UnitIdx);
    end;

    Str(Size:0:1, FormattedNum); { Just one digit after the separator }
    Result := FormattedNum + ' ' + Units[UnitIdx];
end;

{line 359 "phiskyandd.nw"}var
    
{line 878 "phiskyandd.nw"}FirstFileIdx, LastFileIdx : Integer;
MpiRank, MpiSize : Integer;
Idx : Integer;
{line 1206 "phiskyandd.nw"}QuantileArray, OverallQuantileArray : TDoubleArray;
{line 1281 "phiskyandd.nw"}CurFileIdx : Integer;
{line 1560 "phiskyandd.nw"}BinnedMap, HitMap : THealpixMap;
OverallBinnedMap, OverallHitMap : THealpixMap;
{line 1786 "phiskyandd.nw"}FitsFile : TFitsFile;
{line 361 "phiskyandd.nw"}    
{line 692 "phiskyandd.nw"}Configuration : TPhiSkyConfiguration;
{line 742 "phiskyandd.nw"}InputData : TPhiSkyInputData;
{line 975 "phiskyandd.nw"}PhiTod : TPhiSkyTod;

{line 363 "phiskyandd.nw"}begin
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
{line 762 "phiskyandd.nw"}with Configuration do
begin
    Log('Loading the ringsets');
    LoadRingsets(MbRingsetFileNames, InterpolationOrder,
                 InputData.MbRingsets);
    LoadRingsets(SlRingsetFileNames, InterpolationOrder,
                 InputData.SlRingsets);

    Log(Format('Ringsets loaded, %s of memory currently used',
               [MemoryUsed(InputData)]));
end;
{line 1213 "phiskyandd.nw"}SetLength(QuantileArray, (LastFileIdx - FirstFileIdx + 1) * Length(Configuration.Quantiles));
{line 1570 "phiskyandd.nw"}InitMap(Configuration.Nside, Ring, BinnedMap);
InitMap(Configuration.Nside, Ring, HitMap);

InitMap(Configuration.Nside, Ring, OverallBinnedMap);
InitMap(Configuration.Nside, Ring, OverallHitMap);
{line 1248 "phiskyandd.nw"}for CurFileIdx := FirstFileIdx to LastFileIdx do
begin
    with Configuration do
    begin
        try
            ProcessFilePair(PointingFileNames[CurFileIdx],
                            TemperatureFileNames[CurFileIdx],
                            InputData, Configuration.QualityFlagMask,
                            PhiTod, BinnedMap, HitMap);

            AppendQuantiles(PhiTod.PhiSky, PhiTod.Valid, Configuration.Quantiles,
                            QuantileArray, CurFileIdx - FirstFileIdx);

            if Configuration.SaveTods then
                SavePhiTod(ConcatPaths([TodFilePath,
                    ChangeFileExt(ExtractFileName(TemperatureFileNames[CurFileIdx]),
                                  '-phisky.fits')]),
                              PhiTod);
        except
            on E : Exception do
                Log(Format('Unable to process files "%s" and "%s" (%s), skipping…',
                           [PointingFileNames[CurFileIdx],
                            TemperatureFileNames[CurFileIdx],
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
{line 368 "phiskyandd.nw"}        finally
            Mpi.Finalize;
        end;
    except
        on E : Exception do WriteLn('Error: ' + E.message);
    end;
end.
