# ratios 1.2.0

## Major changes

Functions are partly renamed to leave the '.' mostly for generic functions

* function ratio.ds is now called: ratioDT
* function ratio.append_smalles is now called: ratio_append_smallest
* function preparation.DT2 is now called: preparationDT2
* Bug in correction for adhering particles of method 'M1' removed: Calculation of Errors is optional for functions 'Error = TRUE/FALSE'. 

## Minor changes

* at several occasions better testing for correct dimensions of input data sets
* function select.VarsElements has new option 'invert = TRUE/FALSE'.
* function ratioDT has new option 'id.vars'. If given, the function is faster.

# ratios 1.1.1

## Major changes

* function ratio.DT1_DT2 is now called: ratio.ds (for "ratio of data sets").
ratio.DT1_DT2 is kept as alias for backwards compatibility
* the data set UpperCrust is used as default DT2_replace now only in the function CorrectionAdheringParticles, no longer in the function preparation.DT2
* documentation for ratio.ds is improved


# ratios 1.1.0

## Major changes

* function ratio.DT1_DT2: new ratio methods log, ar and cr
* function ratio.DT1_DT2: new option 'vars.ref' for the ratio_type 'ar' and 'alr' where the reference element can be provided. If 'vars.ref' is left empty function still asks interactively for the reference element.
* Description file contains now the info from the ReadMe

## Bug fixes

 * bug in CorrectionAdheringParticles: checking of the classes of the entries is now working correctly
 * selectVarsElements finds element duplicates on each object separately and gives a warning
 * method alr: unnecessary looping replaced by vectored version (much faster now)
 * function preparation.DT2: if 'vars' are provided and vars are element abbreviations function checks now if all 'vars' are in the data sets and if there are more 'vars' than columns it gives a warning
 

# ratios 1.0.0

* Added a `NEWS.md` file to track changes to the package.



