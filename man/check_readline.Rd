% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Subfunctions.R
\name{check_readline}
\alias{check_readline}
\title{check_readline}
\usage{
check_readline(x, myletters)
}
\arguments{
\item{x}{character vector to be checked, e.g. value from \code{\link{readline}}}

\item{myletters}{character vector of entries which should be allowed, e.g. c("yes", "no")}
}
\value{
character vector with one of the values given in \code{myletters}
}
\description{
The function checks if a given character string, e.g. value from \code{\link{readline}}, matches one element of on a given character vector.
Without match it keeps asking for a character string by \code{readline} until 'q' for quit or a matching string is provided.
This is often used to avoid crashing of functions due to wrong input in the options.
}
\details{
Usage is e.g. yesno = check_readline(yesno, c("y", "n")
Now the function will make sure that the variable \emph{x} consists either of "y" or of "n".
}
\examples{
possibleEntries = c("today", "yesterday")
myEntry = "today"
# or try another entry which is different from "today":
# myEntry = readline("Enter any word (without quotes):   ")
y = check_readline(x = myEntry, myletters = possibleEntries)

}
\seealso{
Other sub functions: \code{\link{relError_dataset}},
  \code{\link{select.VarsElements}}
}
\author{
Solveig Pospiech
}
