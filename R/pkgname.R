#' Ratios of environmental data
#'
#' The package provides functions for calculating ratios between two data sets containing environmental data like concentration of elements.
#' Ratios can be calculated either by method "simple", "clr" or "alr".
#' Additionally for plants the amount of adhering particles on plants can be estimated and the element concentration corrected by subtraction.
#'
#' Ratios:
#'
#' Calculating ratios is at the first glance a simple operation but it becomes quickly more complex if the two data sets don't have corresponding rows or columns.
#' A set of functions helps to faster and saver calculate ratios:
#' If for a data set DT1 ratios to a second data set DT2 should be calculated the function \code{\link{preparationDT2}} creates a 'new DT2' with equal number of rows and corresponding columns to DT1 with entries and mean of entries from DT2.
#' Errors are calculated by the function \code{\link{relError_dataset}} as well for each data set as for the ratios, too.
#' The function ratio.DT1_DT2 provides methods for six different types of ratios:
#' 1. simple ratios
#' 2. log ratios
#' 3. ar ratios
#' 4. alr ratios
#' 5. cr ratios
#' 6. clr ratios
#'
#' clr and alr ratios are developed from the clr (center logarithmic transformation) and alr (additive logarithmic transformation) concept introduced by Aitchison in 1986 for compositional data, which are data constrained by a constant sum like concentrations of elements.
#' Hence especially these two methods might be of interest if DT1 and DT2 contain compositional data.
#' This is probably the case for most environmental data.
#' The methods 'ar' and 'cr' are the same as 'alr' and 'clr', but without the logarithm.
#'
#'
#' Correction for Adhering Particles
#'
#' Exact and reproducible analysis of element concentrations in plant tissue is the basis for many research fields such as environmental, health, phytomining, agricultural or provenance studies.
#' Unfortunately plant samples collected in the field will always contain particles on their tissue surfaces such as airborne dust or soil particles.
#' If not removed these particles may induce a bias to the element concentrations measured in plant samples.
#' The influence of adhering particles on element concentration in plants is negligible for elements which have a much higher concentration in the plant tissues compared to the adhering material.
#' This is the case for most main or minor nutrient elements such as P, K, Ca, Mg, S, Mn, B, Mo, Zn or Cu.
#' But elements with typically very low concentrations in plant tissue such as Al, Co, Fe, Li, Ni, Ti, Sc, Zr, or REEs, may show significantly altered concentrations measured in the plants due to adhering particles.
#' Mitchell (1960) proposed that elements with concentration ratios of soil to plant above 100 might show biased concentrations in the measured samples.
#'
#' Reducing the impact of adhering particles on trace element concentration in plants is crucial in order to be able to compare elemental composition of plants, e.g. between sampling periods or slightly different sampling methods or for biomonitoring studies.
#' It is also important in studies for plant nutrition to calculate the real uptake of an element by a plant, e.g. for phytoremediation/phytomining or in studies on the trace elements Co, Ni, Mn and Mo.
#'
#' Based on the model that the analyzed plant material is a mixture of plant tissue and a very minor amount of adhering particles we developed a general methods to calculate a correction term for adhering material.
#'
#' The function \code{\link{Correction.AdheringParticles}} provides three different methods to calculate the influence of adhering particles in order to obtain the element concentrations in plants resulting only from uptake.
#'
#' For further reading and details please refer to the publication:
#' Pospiech, S., Fahlbusch, W., Sauer, B., Pasold, T., & Ruppert, H. (2017). Alteration of trace element concentrations in plants by adhering particles–Methods of correction. Chemosphere, 182, 501-508.
#'
"_PACKAGE"

#' @import data.table
#' @importFrom stringr str_trunc str_trim
#' @importFrom stats mad median
NULL
