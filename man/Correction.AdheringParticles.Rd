% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CorrectionAdheringParticles.R
\name{Correction.AdheringParticles}
\alias{Correction.AdheringParticles}
\title{Correction.AdheringParticles}
\usage{
Correction.AdheringParticles(DT1, DT2 = NULL, vars = NULL,
  vars_ignore = c("As", "Se", "Sn", "V", "Be", "Ge", "Pt"), method, element,
  id.vars, group1.vars, group2.vars, var_subgroup, offset = 0,
  use_only_DT2 = TRUE, DT2_replace = NULL, Errors = TRUE,
  return_as_list = TRUE, negative_values = FALSE,
  set_statistical_0 = FALSE, Error_method = "gauss", STD_DT1 = STD_Plant,
  STD_DT2 = STD_Soil, minNr_DT1 = 100, minNr_DT2 = 100)
}
\arguments{
\item{DT1}{data.frame or data.table, samples in rows and variables in columns}

\item{DT2}{data.frame or data.table, samples in rows and variables in columns.}

\item{vars}{optional, character vector of column names of DT1 and DT2, default is function \code{\link{select.VarsElements}}.
Please make sure the columns given in \code{vars} are of class numeric.}

\item{vars_ignore}{character vector of column names, only for 'method 3'.
These variables are ignored for calculating the median of amount of DT2 (\emph{x}) in 'method 3'.
Please note: the functions returns corrected values for these columns because they are only ignored for calculating the median of \emph{x}.
Default is "As", "Se", "Sn", "V", "Be", "Ge" and "Pt". Please see Details for further explanation.}

\item{method}{characters (no character vector!, please give m3 instead of "m3") denoting the method. Options are \emph{m1}, \emph{m2} and \emph{m3} and \emph{subtr}.
Default is m3. Please see details.}

\item{element}{string, only for method 1. Denotes the column with which amount of DT2 (\emph{x}) is to be calculated.}

\item{id.vars}{column with unique (!) entries for each row. Class can be integer (corresponding row numbers) or character (e.g. sample IDs).
If missing, all columns but \code{vars} will be assigned to it.
Please note: Function is faster and more stable if \code{id.vars} is provided.}

\item{group1.vars}{character vector, column name(s) for subsetting DT1 and DT2}

\item{group2.vars}{optional, column name for subsetting DT1 and DT2 if some entries in \code{group1.vars} are empty.}

\item{var_subgroup}{optional, character vector of one column name of DT1. This option affects the only the error calculation, hence it is ignored if \code{Errors} is set to FALSE.
If provided, DT1 is split into subsets by \code{group1.vars} \emph{and} 'var_subgroup' and the error will calculated for each of these subset.
Please read in the Details for further information.}

\item{offset}{numeric, default is 0. The offset diminishes the subtracted amount of DT2 \emph{x}: x = x - offset.
If used with m2 all concentrations will stay > 0.
Reasonable offset is e.g. offset = 0.0001}

\item{use_only_DT2}{logical, default is FALSE. If there are not enough DT2 data of the location should the DT2s of the region be used? If the \code{use_only_DT2} is set to FALSE then the Upper Crust is used for the correction.}

\item{DT2_replace}{optional, if a DT1 sample does not have DT2 data of the corresponding location with this option you can define which data you would like to use as DT2.
Default is the build-in data set UpperCrust (geochemical composition of the earth's upper crust).
If you would like to have something else, please provide a named vector/ one-row data.table with values used instead of DT2.}

\item{Errors}{logical, should absolute errors get calculated appended to the list - output? Default is FALSE.
If Errors are set to TRUE it overrides the option \code{return_as_list} and always returns a list.}

\item{return_as_list}{logical, should the result get returned as list? Default is FALSE.}

\item{negative_values}{logical, should negative values be returned? If set to FALSE negative values are set to 0. Default is FALSE.}

\item{set_statistical_0}{logical, only for method 3.
Should all values of the variables contributing to the median of \emph{x} be set to 0? Default is FALSE.}

\item{Error_method}{method with which the error should be calculated. At the moment you can choose between "gauss" (default) and "biggest".
See Details for explanation.}

\item{STD_DT1}{optional, data.frame or data.table object for calculating errors for DT1, e.g. the standards. Please see Details. If left empty a default of 5.2\% relative error is used.}

\item{STD_DT2}{optional, data.frame or data.table object for calculating errors for DT2, e.g. the standards. Please see Details. If left empty a default of 5.2\% relative error is used.}

\item{minNr_DT1}{minimum numbers of samples/observations in DT1 for calculating a relative error of observations.
If the number of observations of DT1 is smaller than \code{minNr_DT1} the error is calculated via the data set \code{STD_DT1}.
Default is 50.}

\item{minNr_DT2}{minimum numbers of samples/observations in DT2 for calculating a relative error of observations.
If the number of observations of DT1 is smaller than \code{minNr_DT2} the error is calculated via the data set \code{STD_DT2}.
Default is 50.}
}
\value{
data.frame (or data.table if DT1 is data.table) according to \code{method}.
}
\description{
Suppose element data of one data set (DT1) are biased because the concentrations are the result of a mixture of two substances, of which one substance are the element concentrations of DT2.
In order to correct DT1 to \eqn{DT_{corrected}} a fraction of DT2 has to be subtracted from DT1.
The basic equation for the correction is:
\deqn{
DT_{corrected}=\frac{DT1 - x * DT2}{1 - x}
}
whereof \emph{x} is the amount of DT2 to be subtracted.

The function is written for the case that \emph{x} is unknown.
To calculate \emph{x} the condition is that in \eqn{DT_{corrected}} at least one element concentration is zero or known.
Suppose \eqn{vars_{i}} has a very low concentration, close to zero, in \eqn{DT_{corrected}}: \eqn{DT_{corrected}[vars_{i}]=0}, then:
\deqn{
x =  \frac{DT1[vars_{i}]}{DT2[vars_{i}]}
}

The function was developed for the use to correct plant concentrations for adhering particles:
Exact and reproducible analysis of element concentrations in plant tissue is the basis for many research fields such as environmental, health, phytomining, agricultural or provenance studies.
Unfortunately plant samples collected in the field will always contain particles on their tissue surfaces such as airborne dust or soil particles.
If not removed these particles may induce a bias to the element concentrations measured in plant samples.

For full description of the calculations and the background of correction plants for adhering particles please refer to:

Pospiech, S., Fahlbusch, W., Sauer, B., Pasold, T., & Ruppert, H. (2017). Alteration of trace element concentrations in plants by adhering particles–Methods of correction. \emph{Chemosphere, 182, 501-508}.
and the section Details.
}
\details{
The main option of this function is the \code{method} which determines how the amount of DT2 to be subtracted, the \emph{x}, is going to be calculated.
There are four options:
\itemize{
\item Method 1: calculate \emph{x} via a fixed element
\item Method 2: calculate \emph{x} via the element with the smallest ratio between DT1[vars] and DT2[vars]
\item Method 3: calculate \emph{x} via the median of several, very small ratios between DT1[vars] and DT2[vars]
\item Method subtr: calculate the concentrations for \eqn{x * DT2[vars]}
}

To Method 1:
For example using Ti as \code{element} \eqn{DT_{corrected}} is calculated with \eqn{ x = DT1[Ti]/DT2[Ti]}.
Typical elements for the option \code{element} are e.g. Ti, Al, Zr, Sc, ...
This will eventually lead to negative concentrations for some elements.

To Method 2:
This method subtracts the smallest possible content of DT2 from DT1 (smallest \emph{x}).
For each row/sample the element with the smallest \emph{x} of all ratios \eqn{ x = DT1[vars]/DT2[vars]} of each sample is taken as \emph{element},
hence every sample is corrected based on a different \emph{element}.
With this method there are no negative concentrations.

To Method 3:
In order to reduce the uncertainty of the content of DT2 in DT1 (\emph{x}) based on only one element as in method 1 and 2 an average of several \emph{x} of elements can be calculated.
With \eqn{\Delta x} being the absolute error of \emph{x} the median is calculated by all \emph{x} of elements which values \eqn{ x - \Delta x} are smaller than \eqn{ x_{smallest} + \Delta x_{smallest}}.
The value of the median \eqn{\bar{x}} is then used as \emph{x}.
This will eventually lead to negative concentrations for some elements.
Because statistically the \emph{x} of all elements, which error overlaps the error of the element with smallest \emph{x}, are indistinguishable we suggest to set all elements contributing to \eqn{\bar{x}} to zero, because these small values should not be interpreted:
Set option \code{set_statistical_0} to TRUE.

It is advisable to exclude elements with a huge error margin in the option \code{vars_ignore} because they could severely increase the median \eqn{\bar{x}} by "opening" the window of error-ranges for many elements with significantly higher ratios.
This could lead to an unnatural high median \eqn{\bar{x}} resulting into an overcorrection.

If option \code{id.vars} is provided the functions prints the 'group1.vars' and 'id.vars' of the sample.

For examples and more information please refer to:
Pospiech, S., Fahlbusch, W., Sauer, B., Pasold, T., & Ruppert, H. (2017). Alteration of trace element concentrations in plants by adhering particles–Methods of correction. \emph{Chemosphere, 182, 501-508}.
}
\seealso{
Other ratio functions: \code{\link{preparationDT2}},
  \code{\link{ratioDT}},
  \code{\link{ratio_append_smallest}}
}
\author{
Solveig Pospiech
}
