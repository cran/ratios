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



