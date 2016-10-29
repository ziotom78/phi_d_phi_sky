# Calibration of phi_D and phi_sky

The full documentation for this project is in the PDF file [phiskyandd.pdf](https://github.com/ziotom78/phi_d_phi_sky/blob/master/phiskyandd.pdf).

# Requirements

The code has been written in Pascal. You'll need [Free
Pascal](http://freepascal.org/) to compile and run it. Other requirements:

- Noweb (https://www.cs.tufts.edu/~nr/noweb/)
- LaTeX
- CFITSIO (http://heasarc.gsfc.nasa.gov/fitsio/fitsio.html)
- Make (I used (GNU Make)[https://www.gnu.org/software/make/], but I think it
  might work even with BSD Make)

# Instructions

Run `make` to create the source code from the Noweb files, compile it to two
binary programs (`phid` and `phisky`) and create a PDF containing the full code
with the documentation.

# License

The code is released under a MIT license.

# How to cite this code

If you want to cite this code, please provide a link to the URL
https://github.com/ziotom78/phi_d_phi_sky and cite the following papers:

	@article{ refId0,
		author = {{Planck Collaboration}},
		title = {Planck 2013 results. V. LFI
		  calibration},
		DOI= "10.1051/0004-6361/201321527",
		url= "http://dx.doi.org/10.1051/0004-6361/201321527",
		journal = {A&A},
		year = 2014,
		volume = 571,
		pages = "A5",
	}

	@article{ refId0,
		author = {{Planck Collaboration}},
		title = {Planck 2015 results},
		DOI= "10.1051/0004-6361/201526632",
		url= "http://dx.doi.org/10.1051/0004-6361/201526632",
		journal = {A&A},
		year = 2016,
		volume = 594,
		pages = "A5",
	}

The first one refers to `phid`, the second one to `phisky`.
