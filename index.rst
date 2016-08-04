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

Data
====

In the winter of 2015, Dave Monet provided the LSST Simulations team with a
catalog of bright stars for use in simulating the distribution of guide stars.
The schema of this raw catalog is available in the file
``software/README.sstrc5`` of
`this github repository <https://github.com/lsst-sims/smtn-006>`_.  The catalog
consisted of 645 million stars with RA, Dec, proper motion and magnitudes taken
from one of 10 provenance catalogs as detailed in the README: USNO-B4, 2MASS,
2MASS Extended Source Catalog, WISE, PPMXL, SDSS DR7, PS1 Ubercal, Hipparcos +
Tycho 2, APASS DR7+8, or GCVS and Lepine-Shara.

Fitting to LSST magnitudes
==========================

In order to ingest these stars into the LSSTCATSIM database hosted at the
University of Washington, it was necessary to fit each of these stars to an SED
from the CatSim SED library and assign an E(B-V) value to each.  Because none
of the provided catalogs included parallax information, we decided to fit E(B-V)
values heuristically.  For each SED in ``sims_sed_library``, we calculated
the catalog magnitudes (those contained in the original catalog provided by
Dave Monet) along a grid of E(B-V) values.  We then iterated over the stars in
the input catalog, finding the best-fit combination of [SED, E(B-V) value].
We took this E(B-V) value to be the truest amount of dust extinction between
the Earth and the star.  The prescription for performing this fit is as follows
(all of the software reference below can be found in the ``software/fitting/``
directory of `this github repository <https://github.com/lsst-sims/smtn-006>`_)

1. Run the script ``sed_mag_calc.py``.  This will crawl through all of the
   stellar SEDs in `sims_sed_library` and, for each, will calculate the necessary
   magnitudes along our grid of E(B-V).  This grid of E(B-V) and magnitude values
   will be saved to the file `magnitude_grid.txt`. Note: The E(B-V) grid used is
   specified in the first few lines of the ``if __name__ == "main"`` block of the
   script.  It is currently: 0.001 <= E(B-V) < 0.01 in 0.001 steps;
   0.01 <= E(B-V) < 0.3 in 0.01 steps; 0.3 <= E(B-V) < 7.0 in 0.1 steps.

2. Run the script ``ebv_independent_calc.py``.  This will crawl through all of
   the stellar SEDs in `sims_sed_library` and, for each, will calculate the
   E(B-V)-independent data (normalizing magnitude, Teff, metallicity, log(g) and
   unextincted SDSS magnitudes) required by CatSim.  This data will be written to
   the text file ``ebv_independent_data.txt``.

3. Assemble a text file containing a list of the raw ``.csv.gz``files provided
   by Dave Monet that you would like to fit (this is done so that multiple instances
   of step 5 can be run in parallel without stepping on each other).

4. Compile ``fit_bright_stars.cpp`` like ``g++ -o fit fit_bright_stars.cpp``

5. Run ``fit -m magnitude_grid.txt -e ebv_independent_data -i input_list -o output_dir``.
   This will loop over all of the ``.csv.gz`` files specified in your ``input_list``
   (see step 3) and, for each, produce a text file ``output_dir/[csv name]_ebv_grid_fit.txt``
   containing all of the data needed for ingest into CatSim.

The data columns produced by ``fit_bright_stars.cpp`` are as follows

* htmid -- a place holder for an integer calculated by the database on ingest
* star_id -- the unique integer used to identify the star in the input catalog
  (from input catalog)
* ra -- in degrees (from input catalog)
* dec -- in degrees (from input catalog)
* mura -- RA proper motion in milliarcseconds per year (from input catalog)
* mudec -- Dec proper motion in milliarcseconds per year (from input catalog)
* lon -- galactic longitude in degrees (calculated by `fit_bright_stars.cpp`)
* lat -- galactic latitude in degrees (calculated by `fit_bright_stars.cpp`)
* sed -- the name of the SED in ``sims_sed_library`` associated with the star
  (from the fit)
* magnorm -- the observed magnitude (before dust extinction) of the SED in the
  ImSim bandpass (a delta function at 500 nm; from the fit)
* flux_factor -- the multiplicative factor by which to normalize the SED to get
  observed magnitudes (degenerate with magnorm; from the fit)
* E(B-V) -- best fit value (from the fit)
* Teff -- effective temperature of the star (read off of SED file name)
* [Fe/H] -- metallicity of the star (read off of SED file name)
* log(g) -- surface gravity of the star (read off of SED file name)
* lsst_[]_noatm -- observed LSST magnitudes above the atmosphere (from the fit)
* lsst_[]_atm -- observed LSST magnitudes through standard atmosphere at
  airmass=1.2 (from the fit)
* sdss_[](ext) -- observed SDSS magnitudes with dust extinction applied
  (from the fit)
* sdss_[](raw) -- SDSS magnitudes without dust extinction applied
* color_residual -- the RMS error in the colors used to fit the star to the
  SED (from the fit)
* file_name -- the name of the ``.csv.gz`` file from which this star was read
  (together, file_name and star_id uniquely specify the star; from input catalog)

Schema in CatSim database
=========================

The files produced by ``fit_bright_stars.cpp`` are ingested into the CatSim
database and manipulated to match the schema of the other star tables on the
database.  This database can be accessed by following the instructions `here
<https://confluence.lsstcorp.org/display/SIM/Accessing+the+UW+CATSIM+Database>`_.
The table containing these stars is called `bright_stars`.  The columns contained
in that table are:

* htmid -- an integer which helps with spatial queries on the table

* simobjid -- an integer uniquely identifying the star within this table;
  this is different from the star_id column referenced above.  Together,
  htmid and simobjid uniquely identify each star.

* catlaogid -- the same as star_id above (the integer identifying the star
  in the input catalog)

* source_file -- the name of the ``.csv.gz`` file from which the star was read.
  Together, source_file and catalogid uniquely identify each star.

* ra -- in degrees

* decl -- in degrees (note: it is called decl, not dec)

* mura -- the RA proper motion in milli-arcseconds

* mudecl -- the decl proper motion in milli-arcseconds

* gal_l -- galactic longitude in degrees

* gal_b -- galactic latitude in degrees

* sedfilename -- the name of the associated SED in ``sims_sed_library``

* mag_norm -- the PhoSim normalizing magnitude necessary to match the SED to
  the star's magnitudes

* flux_factor -- the multiplicative factor by which you scale the SED to match
  the star's magnitudes (degenerate with mag_norm

* ebv -- the best-fit E(B-V) value associated with the star and SED

* t -- the effective temperature of the SED

* feh -- the metallicity associated with the SED

* logg -- the log of the surface gravity associated with the SED

* umag, gmag, rmag, imag, zmag, ymag -- apparent LSST magnitudes of the star,
  including Milky Way dust extinction and a standard atmosphere at
  airmass=1.2 (calculated from the SED, E(B-V) pair)

* umag_noatm,... -- apparent LSST magnitudes of the star, including Milky Way dust
  extinction, above the atmosphere.

* newSDSSu, newSDSSg, newSDSSr, newSDSSi, newSDSSz -- the apparent SDSS
  magnitudes of the star, including Milky Way dust extinction (calculated from
  the SED, E(B-V) pair)

* sdssu, sdssg, sdssr, sdssi, sdssz, -- the apparent SDSS magnitudes
  negelecting Milky Way dust extinction

* residual -- the RMS colorresidual between the input stellar colors
  and the colors of the best-fit SED and E(B-V) pair (in the input catalog
  magnitudes measured for that star).

* cx, cy, cz -- the Cartesian coordinates of the star on a hypothetical unit
  sphere (like htmid, these exist to facilitate spatial searching of the table).


Validating the Catalog
======================

Scripts to generate plots useful for validating this catalog are provided in the
``software/validation/`` directory of `this github repository <https://github.com/lsst-sims/smtn-006>`_.  Because the catalog is so large, these
scripts come in groups that work together: some scripts read in the catalog, aggregate
useful quantities, and output those quantities as text files; other scripts read in those
text files and use Matplotlib to produce plots.  The useful groups of scripts are:

* ``dec_vs_mag.py`` queries the CatSim database for all stars between -20.0 < RA < 20.0
  and produces a density plot of mangitude versus Dec in each of the LSST bands.

* ``stellar_density_control_arrays.py reads in the original ``.csv.gz`` files and compiles
  them into HEALPIX maps of number density in 0.5 magnitude bins in the input
  (u, g, r, i, z, y) bands.  These HEALPIX maps are output to text files as simple numpy
  arrays of the number of stars in each HEALPIXel.  ``stellar_density_get_arrays.py``
  queries the CatSim database and assembles the stars into HEALPIX maps of number density
  in 0.5 magnitude bins in each of the LSST bands.  ``stellar_density_comparsions.py`` reads
  in the text files produced by the first two scripts and, for each magnitude bin, plots
  the number density in the input catalog next to the number density in the CatSim database.
  Note: these will be different since, in the input catalog, not every star has every
  magnitude measured.

* ``validate_magnitudes.py`` loops over all of the stars and compiles the number of stars
  in 0.1 magnitudes bins in both (input magnitude, magnitude residual) space as well as
  (magnitude residual, color residual) space.  These grids are output as text files.
  ``plot_magnitude_grids.py`` reads in these text files and produces density plots in both
  of those parameter spaces, as well as cumulative distributions of stars as a function of
  magnitude residual with different cuts applied to color residual.
