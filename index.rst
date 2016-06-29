..
  Content of technical report.

  See http://docs.lsst.codes/en/latest/development/docs/rst_styleguide.html
  for a guide to reStructuredText writing.

  Do not put the title, authors or other metadata in this document;
  those are automatically added.

  Use the following syntax for sections:

  Sections
  ========

  and

  Subsections
  -----------

  and

  Subsubsections
  ^^^^^^^^^^^^^^

  To add images, add the image file (png, svg or jpeg preferred) to the
  _static/ directory. The reST syntax for adding the image is

  .. figure:: /_static/filename.ext
     :name: fig-label
     :target: http://target.link/url

     Caption text.

   Run: ``make html`` and ``open _build/html/index.html`` to preview your work.
   See the README at https://github.com/lsst-sqre/lsst-report-bootstrap or
   this repo's README for more info.

   Feel free to delete this instructional comment.

:tocdepth: 1

In the winter of 2015, Dave Monet provided the LSST Simulations team with a
catalog of bright stars for use in simulating the distribution of guide stars.
The schema of this raw catalog is available in the file
``software/README.sstrc5`` of
`this github repository <https://github.com/lsst-sims/smtn-006>`_.  The catalog
consisted of 645 million stars with RA, Dec, proper motion and magnitudes taken
from one of 10 provenance catalogs as detailed in the README: USNO-B4, 2MASS,
2MASS Extended Source Catalog, WISE, PPMXL, SDSS DR7, PS1 Ubercal, Hipparcos +
Tycho 2, APASS DR7+8, or GCVS and Lepine-Shara.  In order to ingest these stars
into the LSSTCATSIM database hosted at the University of Washington, it was
necessary to fit each of these stars to an SED from the CatSim SED library and
assign an E(B-V) value to each.  Because none of the provided catalogs included
parallax information, we decided to fit E(B-V) values heuristically.  Using the
python script ``software/sed_mag_calc.py`` in `the same github repository
<https://github.com/lsst-sims/smtn-006>`_, we calculated the magnitudes of every
stellar spectrum in ``sims_sed_library`` in every possible bandpass from the 10
catalogs above at values of E(B-V) from 0.001 to 7.0 in 0.1 magnitude
increments (we also included values from 0.001 to 0.009 in 0.001 magnitude
increments to account for objects with very little extinction).  We stored this
grid of magnitudes, SEDs, and E(B-V) values in a text file.  We also used the
script ``software/get_raw_sdss_mags.py`` to calculate the unextincted SDSS
magnitudes of every stellar spectrum in ``sims_sed_library`` (since the CatSim
database schema expects that data as well).  We then used the C++ code
``software/fit_bright_stars.cpp`` to find the best-fitting SED, E(B-V) pair
in N-color space (where N is the maximum number of available colors given the
photometry provided for the star) for each star in the catalog.  We took this
combination to be the true SED and E(B-V) value for the star.  Because the
names of the SED files in ``sims_sed_library`` include the effective
temperatures, metallicities, and surface gravities of the simulate spectra, we
were able to add this data to the table.  Thus, the final schema for the
``bright_stars`` table on the LSSTCATSIM database is as follows:

``bright_stars`` schema
-----------------------
* htmid (an integer which helps with spatial queries on the table)

* simobjid (an integer uniquely identifying the star; taken from the input
  catalog)

* ra in degrees

* decl in degrees

* mura is the RA proper motion in milli-arcseconds

* mudecl is the decl proper motion in milli-arcseconds

* gal_l is galactic longitude in degrees

* gal_b is galactic latitude in degrees

* sedfilename is the name of the associated SED in ``sims_sed_library``

* mag_norm is the PhoSim normalizing magnitude necessary to match the SED to
  the star's magnitudes

* flux_factor is the multiplicative factor by which you scale the SED to match
  the star's magnitudes

* ebv is the best-fit E(B-V) value associated with the star and SED

* t is the effective temperature of the SED

* feh is the metallicity associated with the SED

* logg is the log of the surface gravity associated with the SED

* umag, gmag, rmag, imag, zmag, ymag are apparent LSST magnitudes of the star,
  including Milky Way dust extinction (calculated from the SED, E(B-V) pair)

* newSDSSu, newSDSSg, newSDSSr, newSDSSi, newSDSSz are the apparent SDSS
  magnitudes of the star, including Milky Way dust extinction (calculated from
  the SED, E(B-V) pair)

* sdssu, sdssg, sdssr, sdssi, sdssz, are the apparent SDSS magnitudes
  negelecting Milky Way dust extinction

* residual is the least square error residual between the input stellar colors
  and the colors of the best-fit SED and E(B-V) pair
