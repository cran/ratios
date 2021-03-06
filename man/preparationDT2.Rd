% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Ratios.R
\name{preparationDT2}
\alias{preparationDT2}
\title{preparationDT2}
\usage{
preparationDT2(DT1, DT2, vars = NULL, group1.vars, group2.vars = NULL,
  Errors = FALSE, use_only_DT2 = FALSE, DT2_replace = NULL, minNr = 7,
  STD = NULL, return_as_list = FALSE)
}
\arguments{
\item{DT1}{data.frame or data.table, samples in rows and variables in columns}

\item{DT2}{data.frame or data.table, samples in rows and variables in columns.}

\item{vars}{optional, character vector of column names of DT1 and DT2, default is function \code{\link{select.VarsElements}}.
Please make sure the columns given in \code{vars} are of class numeric.}

\item{group1.vars}{character vector, column name(s) for subsetting DT1 and DT2}

\item{group2.vars}{optional, column name for subsetting DT1 and DT2 if some entries in \code{group1.vars} are empty.}

\item{Errors}{logical, should absolute errors get calculated appended to the list - output? Default is FALSE.
If Errors are set to TRUE it overrides the option \code{return_as_list} and always returns a list.}

\item{use_only_DT2}{logical, default is FALSE. If there are not enough DT2 data of the location should the DT2s of the region be used? If the \code{use_only_DT2} is set to FALSE then the Upper Crust is used for the correction.}

\item{DT2_replace}{mandatory if \code{use_only_DT2} is set to FALSE, serves as substitute for DT2 where DT2 has no corresponding rows to DT1.
A named vector or one-row data.table/ data.frame with the all \code{vars} present.
A column for \code{group1.vars} is not necessary.}

\item{minNr}{minimum numbers of samples/observations for calculating a relative error of observations.
If the number of samples of \code{DT2} is smaller than \code{minNr} the error is calculated via the data set STD.}

\item{STD}{data set for calculating the relative errors if in DT2 there are less rows per group than \code{minNr}.
This replacement data set could for e.g. consist of reference standards with repeated measurement for each standard.}

\item{return_as_list}{logical, should the result get returned as list? Default is FALSE.}
}
\value{
data.frame, data.table or a list, controlled by option \code{return_as_list}.
If \code{Errors} is set to TRUE \code{return_as_list} is ignored and return value is always a list.
The list contains one element if \code{Errors} is set to FALSE and two elements if \code{Errors} is TRUE:
[[1]] is data.table or data.frame of corresponding DT2s, [[2]] data.table or data.frame of absolute errors of corresponding DT2s.
}
\description{
The function creates a data frame 'new DT2' from the variables \code{vars} of the data set DT2 with corresponding rows to the data set DT1,
hence 'new DT2' and DT1 have the same number of rows.
The aim is to generate corresponding rows for two data sets with differing dimensions and even differing number of rows for each group of rows.
For example if for one row i in DT1 there are 3 corresponding rows (j,k,l) in DT2 the function calculates for the 'new DT2' for each variable of \code{vars} an average over the rows j,k and l of DT2,
generating only one row corresponding to the row i in DT1.
If on the other hand for row y in DT2 there are 4 corresponding rows in DT1 the 'new DT2' will contain four times the row y of DT2 matching the four rows of DT1.
The column \code{group1.vars} (and optional \code{group2.vars}) determines which rows of DT1 and DT2 are corresponding, so \code{group1.vars} in DT1 is the look-up table for creating the 'new DT2'.
Generally DT1 and DT2 have to have the columns in common which are given in \code{group1.vars}, \code{group2.vars} and \code{vars}.
}
\details{
The data set 'new DT2' is generated according to following rules:
If there is more than one row in DT2 with the same entry for \code{group1.vars} for each column in \code{vars} an average (mean) of these rows of DT2 is calculated.
After this operation there is only one row for each entry value of \code{group1.vars}.
Each row of this averaged DT2 is replicated n times, with n being the number of rows of the subset of DT1 with the corresponding value in \code{group1.vars}.
If there are values in column \code{group1.vars} in DT1 which are not in DT2 and if option \code{use_only_DT2} is set to \code{TRUE} empty rows are generated.
If option \code{use_only_DT2} is set to \code{FALSE}, data from 'DT2_replace' are taken as substitute for DT2 to fill these empty rows.
The default 'DT2_replace' are element concentrations from the UpperCrust (Rudnick, R. L., & Gao, S. 2003. Composition of the continental crust. \emph{Treatise on geochemistry, 3}, 659.)
}
\seealso{
Other ratio functions: \code{\link{Correction.AdheringParticles}},
  \code{\link{ratioDT}},
  \code{\link{ratio_append_smallest}}
}
\author{
Solveig Pospiech
}
