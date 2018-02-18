#' check_readline
#'
#' @author Solveig Pospiech
#'
#' @description The function checks if a given character string, e.g. value from \code{\link{readline}}, matches one element of on a given character vector.
#' Without match it keeps asking for a character string by \code{readline} until 'q' for quit or a matching string is provided.
#' This is often used to avoid crashing of functions due to wrong input in the options.
#'
#' @param x character vector to be checked, e.g. value from \code{\link{readline}}
#' @param myletters character vector of entries which should be allowed, e.g. c("yes", "no")
#'
#' @details Usage is e.g. yesno = check_readline(yesno, c("y", "n")
#' Now the function will make sure that the variable \emph{x} consists either of "y" or of "n".
#'
#' @examples
#' possibleEntries = c("today", "yesterday")
#' myEntry = "today"
#' # or try another entry which is different from "today":
#' # myEntry = readline("Enter any word (without quotes):   ")
#' y = check_readline(x = myEntry, myletters = possibleEntries)
#'
#' @return character vector with one of the values given in \code{myletters}
#'
#' @family sub functions
#' @export

check_readline<-function(x, myletters)
{
  if (!is.character(myletters)) stop("'myletters' must be a character vector")
  if (!is.character(x)) stop("'x' must be a character vector")
  myletters = c(myletters, "q")
  while(sum(!x %in% myletters)  > 0){
    print.noquote(paste("Please enter", paste(myletters[1:(length(myletters)-1)], collapse = ", "), "(without quotes) or", myletters[length(myletters)], "for 'quit':" ))
    x = readline()
    if(x == "q"){stop("Function had been stopped.")}
  }
  return(x)
}

#' select.VarsElements
#'
#' @author Solveig Pospiech
#'
#' @description The function returns a character vector of element abbreviations if the input object contained variables with element abbreviations.
#' Input may be data.frame, matrix, character vector or named numeric.
#' There are two options to use this functions:
#' \itemize{
#' \item A) only one object
#' \item B) with two objects
#' }
#' For A) the function checks for the pattern of element abbreviations, e.g. Al, S, Ca, etc.
#' For B) the function checks for element abbreviations which are present in both objects. E.g. if x = c("Al", "Ba", "Ca") and y = c("Ba", "K", "Th") the return value will be "Ba".
#' The resulting character vector is without duplicated entries, e.g. x = c("N", "P", "S", "S") results into c("N", "P", "S").
#'
#' @param x data.frame, character vector or named numeric containing element abbreviations as variables
#' @param y optional, data.frame, character vector or named numeric containing element abbreviations as variables
#' @param invert logical. If TRUE return variable names that do not match an element abbreviation pattern
#'
#' @examples
#' x = c("Al", "Ba", "Ca")
#' y = c("Ba", "K", "Th")
#' select.VarsElements(x, y)
#'
#' myvector = c("Al", "Location", "Date", "S", "Ba", "OH")
#' select.VarsElements(myvector)
#' select.VarsElements(myvector, invert = TRUE)
#'
#' @return character vector of element abbreviations
#'
#' @family sub functions
#' @export

select.VarsElements <- function (x, y, invert = FALSE)
{
  if(length(dim(x)) == 2) zw1 = names(x)
  else{
    if(is.character(x)) zw1 = x
    else{
      if(is.numeric(x)) zw1 = names(x) # for named.numeric
      else stop("x must be of class data.frame, character vector or a named numeric.")
    }
  }

  totake = c(grepl("^[A-Z][a-z]$", zw1) | grepl("^[HBCNOFPSKVYIWU]$", zw1)) # look for element patterns
  if(invert) totake = !totake
  Elements = zw1[totake]
  if(sum(duplicated(Elements))>0){
    print.noquote("")
    print.noquote("PLEASE NOTE: In the (first) object there is for some elements more than one entry. Only first column is taken.")
    print.noquote("Following elements are dropped:")
    print.noquote(Elements[duplicated(Elements)])
    Elements = Elements[!duplicated(Elements)]
    print.noquote("")
  }

  if(!missing(y)){
    if(length(dim(y)) == 2) zw2 = names(y)
    else{
      if(is.character(y)) zw2 = y
      else{
        if(is.numeric(y)) zw2 = names(y) # for named.numeric
        else stop("y must be of class data.frame, character vector or a named numeric.")
      }
    }

    totake = c(grepl("^[A-Z][a-z]$", zw2) | grepl("^[HBCNOFPSKVYIWU]$", zw2)) # look for element patterns
    if(invert) totake = !totake
    Elements2 = zw2[totake]
    if(sum(duplicated(Elements2))>0){
      print.noquote("")
      print.noquote("PLEASE NOTE: In the second object there is for some elements more than one entry. Only first column is taken.")
      print.noquote("Following elements are dropped:")
      print.noquote(Elements2[duplicated(Elements2)])
      Elements2 = Elements2[!duplicated(Elements2)]
      print.noquote("")
    }

    Elements = c(Elements, Elements2)
    Elements = Elements[duplicated(Elements)] # take only elements which are in both
  }

  return(sort(Elements))
}


#' relError_dataset
#'
#' @author Solveig Pospiech
#'
#' @description The function calculates for each observation for every variable 'vars' in 'Data' the relative error by median absolute deviation (\code{\link{mad}}) and median (\code{\link{median}}):
#' \deqn{
#' \delta Data[vars_{i}] = \frac{mad(Data[vars_{i}], na.rm = T)}{median(Data[vars_{i}], na.rm = T)}
#' }
#' The observations (e.g. samples) are subset into groups by the column \code{group1.vars}.
#' The relative error is calculated by 'Data' if there are more than 'minNr' entries for each subset of observations.
#' If there are less observations than 'minNr' for a group in 'Data' than the relative error will be calculated by a replacement data set 'STD',
#' e.g. you could use a data set of standard reference samples measured at the same machine as your samples.
#' If you would like to calculate the relative error of all observations in 'Data' set \code{group1.vars} to the column of your sample ID (column with unique entries) and set minNr = 1.
#'
#' @param Data a data.frame or matrix with samples (observations) as rows.
#' @param minNr minimum numbers of samples/observations for calculating a relative error of observations.
#' If the number of samples of \code{Data} is smaller than \code{minNr} the error is calculated via the data set STD.
#' @param STD data set for calculating the relative errors if in \code{Data} there are less rows per group than \code{minNr}.
#' This replacement data set could for e.g. consist of reference standards with repeated measurement for each standard.
#' @param vars optional, character vector of variables of 'Data' for which the error should be calculated.
#' If left empty the function \code{\link{select.VarsElements}} will try to find element abbreviations in the variables of 'Data' and 'STD' if STD is provided.
#' @param group1.vars character vector of variables in 'Data' for splitting 'Data' into subsets. Error will be calculated for each subset.
#' @param group2.vars optional, if a variable name of 'Data' is given here a second splitting by \code{group1.vars} + \code{group2.vars} into subsets is performed.
#' If for grouping by \code{group1.vars} for one subset there are less entries than 'minNr' the function will look up in the second subset
#' if there are enough entries (> minNr) in the group2 corresponding to group1.
#' For example if \code{group1.vars} = "Month" then \code{group2.vars} = "Year" would fill up the gaps if in one month there had been less than 'minNr' observations.
#'
#' @return data.frame or data.table with relative errors for each observation of 'Data'.
#'
#' @family sub functions
#' @export

relError_dataset <- function(Data,
                             vars,
                             group1.vars,
                             group2.vars = NULL,
                             minNr = 7,
                             STD)
{
  UseMethod("relError_dataset", Data)
}

#' @export
relError_dataset.data.table <- function(Data,
                                        vars,
                                        group1.vars,
                                        group2.vars = NULL,
                                        minNr = 7,
                                        STD)
{
  if(missing(group1.vars)) stop("'group1.vars' is missing. Please provide a column name for 'group1.vars'.
                                If you are in doubt choose the column with e.g. your sample IDs and set minNr = 1.")
  if(!is.numeric(minNr)) stop("'minNr' must be numeric.")
  if(length(minNr)>1) stop("minNr must contain only one single number.")
  if(!missing(STD))
    if(!is.data.table(STD)){
      if(!is.data.frame(STD)){
        stop("'STD' must be of class data.frame or data.table.")
      }
      STD_is_dataframe = c()
      STD = data.table(STD)
    }
  print.noquote("")
  print.noquote(paste("--- Calculate relative error by minimum", minNr, "samples ---"))
  print.noquote("")
  if(missing(STD)){
    if(missing(vars)) vars = select.VarsElements(Data)
    if(length(vars) == 0) stop("Please make sure that 'Data' has column names with element abbreviations
                               or provide other column names in 'vars'.")
  }else{
    if(missing(vars)) vars = select.VarsElements(Data, STD)
    if(length(vars) == 0) stop("Please make sure that 'Data' and 'STD' have both columns with element abbreviations
                               or provide column names exisiting both in 'Data' and 'STD' in 'vars'.")
    if(sum(!vars %in% names(STD))  > 0){
      print.noquote("")
      print.noquote("PLEASE NOTE:
                    Please make sure that the object 'STD' contains all variables/columns given in 'vars'
                    or - if you left 'vars' empty - that all columns of 'Data' named with element abbreviations are also present in 'STD'.
                    Missing 'vars' in STD are going to be ignored for relative error calculation.")
      print.noquote("At the moment following columns are having NA for the errors:")
      print.noquote(vars[!vars %in% names(STD)])
      STD[, (vars[!vars %in% names(STD)]) := NA]
      print.noquote("")
    }
  }

  relError = data.table() # create empty table
  ZW = droplevels(Data)[, .N, by = group1.vars]
  for(j in names(ZW)[-length(names(ZW))]){set(ZW, j = j, value = as.character(ZW[[j]]))} # set all columns except last one to character
  if(!is.null(group2.vars)){ # do the same for group2.vars
    ZW2 = droplevels(Data)[, .N, by = c(group1.vars, group2.vars)]
    for(j in c(group1.vars, group2.vars)){set(ZW2, j = j, value = as.character(ZW2[[j]]))}
  }

  if(length(group1.vars) > 1){ # if the group1.vars contains more than one column character
    group1.vars_original = group1.vars
    group1.vars = paste(group1.vars, collapse = "") # be careful, group1.vars gets replaced
    ZW[, (group1.vars) := apply(.SD, 1, function(x) paste(x, collapse = "")), .SDcols = group1.vars_original]
    Data[, (group1.vars) := apply(.SD, 1, function(x) paste(x, collapse = "")), .SDcols = group1.vars_original]
  }

  if(!missing(STD)){
    if(sum(is.na(ZW[[group1.vars]]))  > 0){
      Data_error = data.table(t(STD[, sapply(.SD, function(x) mad(x, na.rm = T)/median(x, na.rm = T)), .SDcols = vars]))
      Data_error[, (group1.vars) := NA]
      relError = rbindlist(list(relError, Data_error), use.names = T, fill = T)
      ZW = ZW[!is.na(ZW[[group1.vars]])]
    }
  }

  for(i in 1:nrow(ZW)){
    myMessage = paste0("For '", ZW[[group1.vars]][i], "'")
    if(ZW$N[i] <= minNr){
      if(!is.null(group2.vars)){
        if(is.na(ZW2[ZW2[[group1.vars]] %in% ZW[[group1.vars]][i]][[group2.vars]])){
          if(missing(STD)){stop("Please provide an entry for STD, because it is missing for calculating the relative errors.")}
          myMessageP = paste(myMessage, "error was calculated from the object STD")
          Data_error = data.table(t(STD[, sapply(.SD, function(x) mad(x, na.rm = T)/median(x, na.rm = T)), .SDcols = vars]))
          Data_error[, (group1.vars) := ZW[[group1.vars]][i]]
        }else{
          if(sum(ZW2[ZW2[[group2.vars]] %in% ZW2[ZW2[[group1.vars]] %in% ZW[[group1.vars]][i]][[group2.vars]]]$N)  > minNr){ # group2.vars has enough measurements
            myMessageP = paste(myMessage, paste0("error was calculated from '", group2.vars, "' because there hadn't been enough measurements in '", paste0(group1.vars, collapse = "' and '"), "'."))
            setkeyv(Data, group2.vars)
            Data_error = data.table(t(Data[ZW2[[group2.vars]][i], sapply(.SD, function(x) mad(x, na.rm = T)/median(x, na.rm = T)), .SDcols = vars]))
            Data_error[, (group1.vars) := ZW[[group1.vars]][i]]
          }else{
            if(missing(STD)){stop("Please provide an entry for STD, because it is missing for calculating the relative errors.")}
            if(sum(!vars %in% names(STD))  > 0){stop("Please make sure that in for calulating the relative error the object 'STD' has the same columns as the object 'Data'.")}
            myMessageP = paste(myMessage, "error was calculated from the object STD")
            Data_error = data.table(t(STD[, sapply(.SD, function(x) mad(x, na.rm = T)/median(x, na.rm = T)), .SDcols = vars]))
            Data_error[, (group1.vars) := ZW[[group1.vars]][i]]
          }
        }
      }else{
        if(missing(STD)){stop("Please provide an entry for STD, because it is missing for calculating the relative errors.")}
        myMessageP = paste(myMessage, "error was calculated from the object STD")
        Data_error = data.table(t(STD[, sapply(.SD, function(x) mad(x, na.rm = T)/median(x, na.rm = T)), .SDcols = vars]))
        Data_error[, (group1.vars) := ZW[[group1.vars]][i]]
      }
    }else{
      myMessageP = paste0("Error calculated from '", paste0(group1.vars, collapse = "' and '"), "'.")
      setkeyv(Data, group1.vars)
      Data_error = data.table(t(Data[ZW[[group1.vars]][i], sapply(.SD, function(x) mad(x, na.rm = T)/median(x, na.rm = T)), .SDcols = vars]))
      Data_error[, (group1.vars) := ZW[[group1.vars]][i]]
    }
    relError = rbindlist(list(relError, Data_error), use.names = T, fill = T)
    print.noquote(myMessageP)
  }
  setkeyv(relError, group1.vars)
  relError = relError[as.character(Data[[group1.vars]])]
  if(exists("group1.vars_original")){Data[, (group1.vars) := NULL]} # remove the additional column from Data
  if(exists("STD_is_dataframe")){STD = data.frame(STD)}
  return(relError)
}

#' @export
relError_dataset.data.frame <- function(Data,
                                        vars,
                                        group1.vars,
                                        group2.vars = NULL,
                                        minNr = 7,
                                        STD)
{
  Data = data.table(Data)
  relError = relError_dataset.data.table(Data = Data,
                                         minNr = minNr,
                                         vars = vars,
                                         group1.vars = group1.vars,
                                         group2.vars = group2.vars,
                                         STD = STD)
  Data = data.frame(Data)
  relError = data.frame(relError)
  return(relError)
}

#' @export
relError_dataset.default <- function(Data,
                                     vars,
                                     group1.vars,
                                     group2.vars = NULL,
                                     minNr = 7,
                                     STD)
{
  if(missing(Data)){stop("Please provide a dataset of which the error should be calculated")}
  stop("'Data' must be a data.frame or a data.table object.")
}

