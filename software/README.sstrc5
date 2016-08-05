SSTRC5  --  Revision of January 2015

I)  INTRODUCTION

   This revision should give the same astrometry as previous but should
   offer substantial increases in the photometric accuracy.  As with
   SSTRC, it is derived from a merge and extraction from catalogs
   compiled at or available to USNO Flagstaff personnel.  In order of
   compilation (which is the inverse of precedence), the catalogs are:

   USNO-B4:
      B4 is substantially different than USNO-B1.0 in that (a) the
      systematic errors between B1 and 2MASS have been removed, and
      (b) the photographic magnitudes (O,E,J,F,N) have been replaced
      by PS1-based (g,r,i).  The color terms in these new magnitudes
      have been ignored due to the much larger errors in photographic
      photometry.

   2MASS:
      2MASS coordinates replace B4 values and its (J,H,K) photometry is
      appended.

   2MASS Extended Source Catalog:
      Added for completeness.  In the culling process (below), most of
      these entries are removed.

   WISE:
      Based on AllWISE table of fluxes associated with a 2MASS object.
      No WISE fluxes for non-2MASS objects.

   PPMXL:
      Update previous coordinates and motions.

   SDSS DR7:
      Coordinates are used if objects are not it 2MASS, and its values
      for (u,g,r,i,z) are appended.

   PS1 Ubercal (PV2 of Dec 2014):
      Coordinates are used if objects are not in 2MASS or SDSS.  PS1
      (g,r,i,z) replace SDSS values if more colors are available from
      PS1 than from SDSS.

   UCAC4:
      UCAC4 coordinates are used only if the object
      was not in 2MASS, SDSS, or PS1 because the systematics of the
      error distribution are not yet known.  The proper motions replace
      those found earlier, and the (r) magnitude is used only if no other
      (r) magnitude is available.

   Hipparcos + Tycho2:
      Coordinates are used in preference to any other catalog, and its
      (B,V) magnitudes are appended.

   APASS DR7+8:
      Values for (B,V,g,r,i) are used.

   GCVS and Lepine-Shara
      Flags are set if object appears in General Catalog of Variable Stars
      dated Apr 2012 and the full sky version of the Lepine-Shara catalog
      of large proper motion stars.

   The master catalog (PNOMAD) contains all objects from all catalogs.
   When possible, objects are merged on the basis of their (RA,Dec).
   After the master catalog is compiled, a separate culling program
   extracts objects for a specific catalog tailored for a specific
   telescope or sensor system.  A discussion of the culling process
   for SSTRC5 is given below.

   The culled catalog is contained in 1800 files containing CSV lines.
   As requested, NULL entries appear as the string "-99".  The schema
   for the lines is presented below.  Because NULL entries can appear
   only for magnitudes, there is no chance for confusion.

II)  SCHEMA

   The catalog is contained in 1800 files, each labeled for a 0.1 degree
   wide stripe of South Polar Distance (Declination + 90 degrees; 0 to 1799).
   In each file, the objects are sorted in RA.  Each line contains 22
   fields, and the format is CSV.  See Note 5 regarding NULL values.

   Column   Name      Type    Description
   ---------------------------------------------------------------------------
      1     starID   (long)   Unique star number.  Note 1.
      2     ra      (double)  J2000 Right Ascension in degrees.
      3     decl    (double)  J2000 Declination in degrees.
      4     muRA    (float)   Proper motion in RA in arcsec/year.  Note 2.
      5     muDecl  (float)   Proper motion in Dec in arcsec/year.  Note 2.
      6     magB    (float)   Magnitude in B band.
      7     magV    (float)   Magnitude in V band.
      8     magU    (float)   Magnitude in u band.
      9     magG    (float)   Magnitude in g band.
     10     magR    (float)   Magnitude in r band.
     11     magI    (float)   Magnitude in i band.
     12     magZ    (float)   Magnitude in z band.
     13     magY    (float)   Magnitude in y band.
     14     magJ    (float)   Magnitude in J band.
     15     magH    (float)   Magnitude in H band.
     16     magK    (float)   Magnitude in K band.
     17     w1      (float)   Magnitude in WISE W1 band.
     18     w2      (float)   Magnitude in WISE W2 band.
     19     w3      (float)   Magnitude in WISE W3 band.
     20     w4      (float)   Magnitude in WISE W4 band.
     21     magSST  (float)   Magnitude in SST band.  Note 3.
     22     flag    (int)     Flag bits.  Note 4.

   Notes to Schema:

   1) The starID comes from one of the component catalogs.  When examined
      in hexadecimal, the following base values indicate which catalog
      provided the ID.

      0x100000000000    SDSS star
      0x200000000000    SDSS galaxy
      0x300000000000    B4
      0x400000000000    2MASS
      0x500000000000    UCAC
      0x600000000000    Tycho2 or Hipparcos
      0x700000000000    APASS
      0xa00000000000    PS1 Ubercal

   2) Proper motions are in the tangent plane and are given in arcsec/year.
      Conversion to spherical coordinates requires dividing the RA proper
      motion by cos(declination).

   3) magSST is filled with the best available r magnitude.

   4) The flag bits are set during the compilation of the catalog.  When
      examined as a hexadecimal number, the bit values are:

      Value        Description
     ------------------------------------------------------------------------
     0x00000001    Object is in B4 catalog.
     0x00000002    Object is in 2MASS catalog.
     0x00000004    Object is in SDSS DR7 catalog.
     0x00000008    Object is in PPMXL catalog.
     0x00000010    Object is in UCAC4 catalog.
     0x00000020    Object is in Tycho2 catalog.
     0x00000040    Object is in Hipparcos catalog.
     0x00000080    Object is in AllWISE catalog.
     0x00000100    Object is in Lepine-Shara catalog.
     0x00000200    Object is in APASS catalog.
     0x00000400    Object is in PS1 catalog.
     0x00000800    Object is in 2MASS Extended Source catalog.
     0x00100000    CLIP - cull algorithm sets/removes based on this bit.
     0x00200000    PPMXL has large residuals compared to B4
     0x00400000    WORRY - like CLIP but not fatal; something is curious
     0x00800000    g,r,i come from "accurate" photometry catalog.
     0x01000000    Object marked DIRTY in 2MASS catalog flags.
     0x02000000    Object marked SPIKE in B4 catalog flags.
     0x04000000    Object marked CONFUSED in Tycho2 catalog flags.
     0x08000000    B,V come from "accurate" photometry catalog.
     0x10000000    Object passes most sanity tests and is probably useful.
     0x20000000    Object has no r band magnitude so magSST is not defined.
     0x40000000    Object marked as galaxy in SDSS or PS1 catalog.
     0x80000000    Object is in GSVS or Lepine-Shara.

   5) The fields
         starID, ra, decl, muRA, muDecl, flag
      cannot contain NULL values.  In the case of the proper motions,
      NULL values have been set to zero.  Only the magnitude columns can
      contain null values.  As requested, a null value appears at the
      string "-99".

   6) The value in 0x0000f000 is the Astrometry Source catalog, and the
      value in 0x000f0000 is the Photometry Source catalog.  Whereas many
      bits can be set in the 0x00000fff fields for any or all catalogs
      containing this object, the ASRC and PSRC contain the catalog
      whose values are quoted.

      Value    Source
      ---------------
         0      B4 photometry from uncalibrated plates, B4 astrometry
         1      B4 photometry from recalibrated plates, B4 astrometry
         2      2MASS
         3      SDSS
         4      PPMXL
         5      UCAC4
         6      Tycho2
         7      Hipparcos
         8      WISE
         9      Lepine-Shara
        10      APASS
        11      PS1
        12      2MASS Extended Source
     13-15      (not used)

III)  CULLING ALGORITHM

   As noted above, the master catalog contains every object in each of the
   component catalogs.  Once compiled, it should remain reasonably constant
   until new versions of the component catalogs (or new component catalogs)
   become available.  As of this document, it contains 3317407385 lines.

   The culled catalog applies a rule set to the master catalog, and extracts
   the desired subset.  For SSTRC5, the culling algorithm is the following.

   1) Reject if any bad bit is set (CLIP, DIRTY, CONFUSE, SPIKE, NOMAG,
         GALAXY).

   2) Reject if proper motion is >=0.2 arcsec/year.

   3) Require keep if star is in Tycho2 or Hipparcos or APASS.  This is
         an attempt to improve the completeness at the bright end.

   4) Reject is star appears in only one catalog that is not believed to
         be trustworthy.

   5) Cull based on r-band magnitude if defined, or V-band if not.  Reject
         if neither is available.
