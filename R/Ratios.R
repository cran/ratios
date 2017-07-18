#' preparation.DT2
#'
#' @author Solveig Pospiech
#'
#' @description The function creates a data set 'new DT2' from the variables \code{vars} of DT2, which has equal number of rows to DT1.
#' The column \code{group1.vars} (and optional \code{group2.vars}) in DT1 is the look-up table for creating the 'new DT2'.
#' Generally DT1 and DT2 have to have the columns in common which are given in \code{group1.vars}, optional \code{group2.vars} and \code{vars}.
#'
#' The data set 'new DT2' is generated according to following rules:
#' If there is more than one row in DT2 with the same entry for \code{group1.vars} for each column in \code{vars} an average (mean) of these rows of DT2 is calculated.
#' After this operation there is only one row for each entry value of \code{group1.vars}.
#' Each row of this averaged DT2 is replicated n times, with n being the number of rows of the subset of DT1 with the corresponding value in \code{group1.vars}.
#' If there are values in column \code{group1.vars} in DT1 which are not in DT2 and if option \code{use_only_DT2} is set to \code{TRUE} empty rows are generated.
#' If option \code{use_only_DT2} is set to \code{FALSE}, data from 'DT2_replace' are taken as substitute for DT2 to fill these empty rows.
#' The default 'DT2_replace' are element concentrations from the UpperCrust (Rudnick, R. L., & Gao, S. 2003. Composition of the continental crust. \emph{Treatise on geochemistry, 3}, 659.)
#'
#' @param DT1 data.frame or data.table, samples in rows and elements and other information in columns
#' @param DT2 data.frame or data.table of DT2 or crust data, samples in rows and elements and other information in columns.
#' @param DT2_replace optional, if a DT1 sample does not have DT2 data of the corresponding location with this option you can define which data you would like to use as DT2.
#' Default is the UpperCrust. If you would like to have something else, please provide a named vector/ one-row data.table with values used instead of DT2
#' @param group1.vars character vector, column name(s) for subsetting DT1 and DT2
#' @param group2.vars optional, column name for subsetting DT1 and DT2 if some entries in \code{group1.vars} are empty.
#' @param vars optional, character vector of column names of DT1 and DT2, default is function \code{\link{select.VarsElements}}.
#' Please make sure the columns given in \code{vars} are of class numeric.
#' @param use_only_DT2 logical, default is FALSE. If there are not enough DT2 data of the location should the DT2s of the region be used? If the \code{use_only_DT2} is set to FALSE then the Upper Crust is used for the correction.
#' @param Errors logical, should absolute errors get calculated appended to the list - output? Default is FALSE.
#' If Errors are set to TRUE it overrides the option \code{return_as_list} and always returns a list.
#' @param return_as_list logical, should the result get returned as list? Default is FALSE.
#' @inheritParams relError_dataset
#'
#' @return data.frame, data.table or a list, controlled by option \code{return_as_list}.
#' If \code{Errors} is set to TRUE \code{return_as_list} is ignored and return value is always a list.
#' The list contains one element if \code{Errors} is set to FALSE and two elements if \code{Errors} is TRUE:
#' [[1]] is data.table or data.frame of corresponding DT2s, [[2]] data.table or data.frame of absolute errors of corresponding DT2s.
#'
#' @family ratio functions
#' @export

preparation.DT2 <- function(DT1,
                            DT2,
                            vars = NULL,
                            group1.vars,
                            group2.vars = NULL,
                            Errors = FALSE,
                            use_only_DT2 = FALSE,
                            DT2_replace = NULL,
                            minNr = 7,
                            STD = NULL,
                            return_as_list = FALSE)
{
  UseMethod("preparation.DT2", DT2)
}

#' @export
preparation.DT2.data.table <- function(DT1,
                                       DT2,
                                       vars = NULL,
                                       group1.vars,
                                       group2.vars = NULL,
                                       Errors = FALSE,
                                       use_only_DT2 = FALSE,
                                       DT2_replace = NULL,
                                       minNr = 7,
                                       STD = NULL,
                                       return_as_list = FALSE)
{
  print.noquote("--- Preparation DT2 ---")
  print.noquote("")
  if(Errors)
    if(is.null(STD)){stop("Please provide a data set for STD or set Errors = F")}
  if(missing(DT1)) stop("Please provide a first data set, to which the ratio should be calulated.")
  else
    if(!is.data.table(DT1)){
      if(!is.data.frame(DT1)) stop("DT1 must be data.frame or data.table")
      DT1 = data.table(DT1)
    }
  if(missing(group1.vars)){stop("Please provide 'group1.vars'.")}
  if(!is.character(group1.vars)) stop("'group1.vars' must be a character vector.")
  if(!group1.vars %in% names(DT1)) stop(paste("In the data set of the DT1 the column", group1.vars, "is missing"))
  if(!group1.vars %in% names(DT2)) stop(paste("In the data set of the DT2 the column", group1.vars, "is missing"))
  if(!is.null(group2.vars)){
    if(!is.character(group2.vars)) stop("'group2.vars must be a character vector.")
    if(!group2.vars %in% names(DT1)) stop(paste("In the data set of the DT1 the column", group2.vars, "is missing"))
    if(!group2.vars %in% names(DT2)) stop(paste("In the data set of the DT2 the column", group2.vars, "is missing"))
  }
  # check if there are congruenting entries in column 'group1.vars' in DT1 and DT2:
  if(sum(levels(as.factor(droplevels(DT1)[[group1.vars]])) %in% levels(as.factor(droplevels(DT2)[[group1.vars]]))) == 0){
    print.noquote(paste("In the column", group1.vars, "DT1 and DT2 have no entries in common."))
    print.noquote("Hence it is impossible to know which observations from DT2 belong to which observations in DT1.")
    if(use_only_DT2) stop("Please provide as 'group1.vars' a column which exist in both, DT1 and DT2, AND has entries in common.")
    else
      print.noquote("WARNING! DT2 is ignored and only DT2_replace is used!")
      print.noquote("")
  }

  if(is.null(vars)){
    vars = select.VarsElements(DT1, DT2)
    if(length(vars) == 0) stop("please provide for 'vars' a character vector with column names which exist in both DT1 and DT2")
  }else{
    vars_DT1 = select.VarsElements(vars, DT1)
    if(length(vars) != length(vars_DT1)){
      print.noquote("WARNING: Following 'vars' are not in DT1 and going to be ignored:")
      print.noquote(vars[!vars %in% vars_DT1])
      print.noquote("")
      vars = vars[vars %in% vars_DT1]
    }
    vars_DT2 = select.VarsElements(vars, DT2)
    if(length(vars) != length(vars_DT2)){
      print.noquote("WARNING: Following 'vars' are not in DT2 and going to be ignored:")
      print.noquote(vars[!vars %in% vars_DT2])
      print.noquote("")
      vars = vars[vars %in% vars_DT2]
    }
  }
  if(length(vars) == 0) stop("please provide for 'vars' a character vector with column names which exist in both DT1 and DT2")

  isitnumeric <- DT2[,sapply(.SD, class), .SDcols=vars] == "numeric"
  if(all(!isitnumeric)) stop("Please provide at least one variable in'vars' which is numeric.")
  if(!all(isitnumeric)){
    print.noquote("There some columns specified in'vars' which are not numeric:")
    print.noquote(vars[!isitnumeric])
    print.noquote("There are going to be ignored for the following calculations.")
    vars = vars[isitnumeric]
  }
  colsinfo = names(DT1)[!names(DT1) %in% vars]
  # make row-numbering for reordering at the end:
  DT1[, "Rows" := c(1:nrow(DT1))]

  if(!use_only_DT2){
    if(is.null(DT2_replace)){
      DT2_replace = UpperCrust[, vars[vars %in% names(UpperCrust)], with = F]
      if(ncol(DT2_replace) == 0) stop("Given columns names are not columns names of the UpperCrust. Please provide another DT2_replace or different column names.")
    }else{
      if(!is.data.frame(DT2_replace)) stop("'DT2_replace' must be data.frame or data.table")
      if(!is.data.table(DT2_replace)) DT2_replace = data.table(DT2_replace)
      if(all(!vars %in% names(DT2_replace))) stop("Given columns names are not columns names of the DT2_replace. Please provide another DT2_replace or different column names.")
    }
    # add empty columns:
    DT2_replace[, (vars[!vars %in% names(DT2_replace)]) := NA]
  }

  setkeyv(DT1, group1.vars)

  Info2 = DT2[, (names(DT2)[names(DT2) %in% colsinfo]), with = F] # not used further, because it is not so easy to be implemented.
  # The idea is to give the information which are the same for all samples contributing the average back.
  DT2_mean_all = DT2[, lapply(.SD, mean, na.rm = T), .SDcols = vars, by = c(group1.vars, group2.vars)] # means of all DT2s
  DT2_mean_all = DT2_mean_all[rowSums(is.na(DT2_mean_all))!= (length(vars)+2)]

  # all group1.varss, which are not in the DT2s: (that ist important if a second group2.vars is given)
  Variables_not_in_DT2 = levels(as.factor(DT1[[group1.vars]]))[!levels(as.factor(DT1[[group1.vars]])) %in% levels(as.factor(DT2_mean_all[[group1.vars]]))]
  if(length(Variables_not_in_DT2)>0){
    print.noquote("There are entries in DT1 which are not in DT2:")
    print.noquote("")
    DT2Mean_all = data.table()
    for(Variable in Variables_not_in_DT2){
      myMessage = paste(group1.vars, Variable)
      if(use_only_DT2){
        if(is.null(group2.vars)){
          myMessageP = paste0(myMessage, paste0(": For '", Variable, "' no entry in DT2 could be found. Everything is replaced by NA"))
          DT2Mean = data.table()
        }else{
          Variable2 = DT1[Variable, mean(get(vars[1])), by = c(group1.vars, group2.vars)][[group2.vars]] # get the Variable2 of the Variable, vars[1] ist nur dummy für irgendein Element
          if(length(Variable2)>1){
            Variable2 = NA_character_
            myMessageP = paste0(myMessage, paste0(": There has been more than one entry in '", group2.vars, "' for '", group1.vars, "'. No mean caculated"))
          }
          setkeyv(DT2_mean_all, group2.vars)
          if(nrow(DT2_mean_all[Variable2, nomatch = 0]) == 0){ # Variable2 not in DT2s
            if(use_only_DT2){
              myMessageP = paste0(myMessage, paste0(": There has been no entry in '", group2.vars, "' in DT2. No mean caculated"))
              DT2Mean = data.table()
            }else{
              myMessageP = paste0(myMessage, paste0(": There has been no entry in '", group2.vars, "' in DT2. DT2_replace used"))
              DT2Mean = copy(DT2_replace)
              DT2Mean[, (group1.vars) := Variable]
            }
          }else{
            myMessageP = paste0(myMessage, paste0(": Mean of '", Variable2, "' taken as DT2 for '", Variable, "'."))
            setkeyv(DT2, group2.vars)
            DT2Mean = DT2[Variable2, lapply(.SD, mean, na.rm = T), .SDcols = vars]
            DT2Mean[, (group1.vars) := Variable]
          }
        }
      }else{ # make DT2 also from DT2_replace
        myMessageP = paste0(myMessage, c(": DT2_replace taken as DT2"))
        DT2Mean = copy(DT2_replace)
        DT2Mean[, (group1.vars) := Variable]
      }
      print.noquote(myMessageP)
      DT2Mean_all = rbindlist(list(DT2Mean_all, DT2Mean), fill = T, use.names = T)
    }
    DT2_mean_all = rbindlist(list(DT2_mean_all, DT2Mean_all), use.names = T, fill = T)
  }else{
    DT2Mean_all = data.table()
  }
  setkeyv(DT2, group1.vars)
  setkeyv(DT1, group1.vars)

  # if there group1.vars multiple times in group2.vars:
  if(sum(DT2_mean_all[, .N, by = group1.vars]$N>1)>0){
    print.noquote(paste0("Warning: For following entries in '", group1.vars, "' there is more than one combination with '", group2.vars, "' in DT2."))
    if(!is.null(group2.vars)){
      print.noquote(paste0("Only the first entry is taken. Please correct column '", group2.vars, "' in DT2 raw-data and rerun function, if you care for it."))
    }
    print.noquote("Following entries in DT2 are affected:")
    print.noquote(DT2_mean_all[, .N, by = group1.vars][DT2_mean_all[, .N, by = group1.vars]$N>1])
    rowremove = which(duplicated(DT2_mean_all[[group1.vars]]))
    DT2_mean_all = DT2_mean_all[-rowremove, ]
    setkeyv(DT2_mean_all, group1.vars)
  }

  setkeyv(DT2_mean_all, group1.vars)
  DT2_mean_all = DT2_mean_all[as.character(DT1[[group1.vars]])]
  if(sum(is.na(DT1[[group1.vars]]))>0){
    if(use_only_DT2){
      print.noquote("BE AWARE:")
      print.noquote(paste0("There are ", sum(is.na(DT1[[group1.vars]])), " NAs in column '", group1.vars, "' in DT1. These are omitted for ratios."))
      print.noquote("Following entries are empty:")
      print.noquote(DT1[is.na(DT1[[group1.vars]]), .N, by = group1.vars])
      print.noquote("")
    }else{
      # make DT2 also from DT2_replace
      print.noquote(paste0("For missing values in column '", group1.vars, "' DT2_replace taken as DT2"))
      print.noquote("")
      for(j in vars){set(DT2_mean_all, j = j, i = which(is.na(DT2_mean_all[[group1.vars]])), value = DT2_replace[[j]])}
    }
  }

  if(Errors){
    suppressWarnings(STD[, (vars[!vars %in% names(STD)]) := NA])
    print.noquote("")
    print.noquote("Errors for DT2:")
    error_DT2_rel = relError_dataset(Data = DT2_mean_all,
                                 group1.vars = group1.vars,
                                 group2.vars = group2.vars,
                                 STD = STD, minNr = minNr)
    error_DT2 = copy(DT2_mean_all)
    for(kol in vars){set(error_DT2, j = kol, value = error_DT2[[kol]]*error_DT2_rel[[kol]])}
    setkeyv(error_DT2, group1.vars)
  }

  # set roworders:
  DT2_mean_all[, "Rows" := DT1[["Rows"]]]
  setkeyv(DT2_mean_all, "Rows")
  DT2_mean_all[, "Rows" := NULL]
  if(Errors){
    error_DT2[, "Rows" := DT1[["Rows"]]]
    setkeyv(error_DT2, "Rows")
    error_DT2[, "Rows" := NULL]
  }
  setkeyv(DT1, "Rows")
  DT1[, "Rows" := NULL]

  if(Errors){
    return(list("DT2" = DT2_mean_all, "error" = error_DT2))
  }else{
    if(return_as_list)
      return(list("DT2" = DT2_mean_all))
    else
      DT2_mean_all
  }
}

#' @export
preparation.DT2.data.frame <- function(DT1,
                                       DT2,
                                       vars = NULL,
                                       group1.vars,
                                       group2.vars = NULL,
                                       Errors = FALSE,
                                       use_only_DT2 = FALSE,
                                       DT2_replace = NULL,
                                       minNr = 7,
                                       STD = NULL,
                                       return_as_list = FALSE)
{
  DT2 = data.table(DT2)
  myReturn=preparation.DT2.data.table(DT1 = DT1,
                                      DT2 = DT2,
                                      DT2_replace = DT2_replace,
                                      group1.vars = group1.vars,
                                      group2.vars = group2.vars,
                                      vars = vars,
                                      use_only_DT2 = use_only_DT2,
                                      Errors = Errors,
                                      STD = STD,
                                      minNr = minNr,
                                      return_as_list = return_as_list)
  if(return_as_list){
    myReturn$DT2 = data.frame(myReturn$DT2)
    if(Errors)
      myReturn$error = data.frame(myReturn$error)
    return(myReturn)
  }else{
    return(data.frame(myReturn))
  }
}

#' @export
preparation.DT2.default <- function(DT1,
                                    DT2,
                                    vars = NULL,
                                    group1.vars,
                                    group2.vars = NULL,
                                    Errors = FALSE,
                                    use_only_DT2 = FALSE,
                                    DT2_replace = NULL,
                                    minNr = 7,
                                    STD = NULL,
                                    return_as_list = FALSE)
{
  stop("'DT2' must be a data.frame or a data.table object.")
}


#' ratio.DT1_DT2
#'
#' @author Solveig Pospiech
#'
#' @description The function calculates ratio between DT1 and DT2 for all variables specified in \code{vars} by the columns \code{group1.vars} (and optional \code{group2.vars}).
#' Generally DT1 and DT2 have to have the columns in common which are given in \code{group1.vars}, optional \code{group2.vars} and \code{vars}.
#' If DT2 has different number of rows than DT1 a 'new DT2' with corresponding dimensions is prepared by the function \code{\link{preparation.DT2}}.
#' At the moment there are three different options for calculating the ratios:
#' \itemize{
#'     \item "simple"
#'     \item "normalized"
#'     \item "cr"
#'     \item "clr"
#'     \item "alr"
#' }
#' For more details please refer to \code{\link{preparation.DT2}} and section Details.
#'
#' @inheritParams preparation.DT2
#' @param ratio_type character vector of "simple", "log", "ar", "alr", "cr" and "clr".
#' Please refer to details for explanations.
#' @param vars.ref reference variable, one out of \code{vars}. Only for \code{ratio_type} "ar" or "alr".
#' @param STD_DT1 optional, data.frame or data.table object for calculating errors for DT1, e.g. the standards. Please see Details. If left empty a default of 5.2\% relative error is used.
#' @param STD_DT2 optional, data.frame or data.table object for calculating errors for DT2, e.g. the standards. Please see Details. If left empty a default of 5.2\% relative error is used.
#' @param return_all logical, should \emph{all} used data sets be returned as a list? Default is FALSE.
#' If set to TRUE the list contains DT1, DT2, vars, ratios, and optional additional ratios_error, DT1_error and DT2_error.
#' @param return_as_list logical, should the result get returned as list? Default is FALSE.
#' If set to FALSE and \code{Errors} is set to TRUE a column \code{type_of_data} is appended.
#' This option is ignored if option 'return_all' is set to TRUE.
#' @param minNr_DT1 minimum numbers of samples/observations in DT1 for calculating a relative error of observations.
#' If the number of observations of DT1 is smaller than \code{minNr_DT1} the error is calculated via the data set \code{STD_DT1}.
#' Default is 50.
#' @param minNr_DT2 minimum numbers of samples/observations in DT2 for calculating a relative error of observations.
#' If the number of observations of DT1 is smaller than \code{minNr_DT2} the error is calculated via the data set \code{STD_DT2}.
#' Default is 50.
#' @param Error_method method with which the error should be calculated. At the moment you can choose between "gauss" (default) and "biggest".
#' See Details for explanation.
#' @param var_subgroup optional, character vector of one column name of DT1. This option affects the only the error calculation, hence it is ignored if \code{Errors} is set to FALSE.
#' If provided, DT1 is split into subsets by \code{group1.vars} \emph{and} 'var_subgroup' and the error will calculated for each of these subset.
#' Please read in the Details for further information.
#'
#' @details To calculate the ratios the functions internally calls \code{\link{preparation.DT2}} to create a data set 'new DT2' from the variables \code{vars} of DT2, which has equal number of rows to DT1.
#' Then the division is done by the now corresponding data sets by the method given in 'ratio_type'.
#'
#' The method "simple" is a simple division between DT1 and DT2:
#' \deqn{
#' \frac{DT1[vars]}{DT2[vars]}
#' }
#' The method "log" is the logarithm of the simple ratio:
#' \deqn{ ln \left( \frac{DT1[vars]}{DT2[vars]} \right)}
#'
#' The methods "ar" and "alr" normalize all ratios to one reference column:
#' ar:
#' \deqn{
#' \frac{DT1[vars_{i}]}{DT2[vars_{i}]} * \frac{DT2[vars_n]}{DT1[vars_n]}_{i=1,\dots, n, \dots, D}
#' }
#' alr:
#' \deqn{
#' ln \left(\frac{DT1[vars_{i}]}{DT2[vars_{i}]} * \frac{DT2[vars_n]}{DT1[vars_n]}\right)_{i=1,\dots, n, \dots, D}
#' }
#'
#' The methods "cr" and "clr" normalize all ratios to the geometric mean of all columns included by \code{vars}:
#' "cr" is calculated by:
#' \deqn{
#' \frac{DT1[vars_{i}]}{DT2[vars_{i}]} * \frac{g(x)^{DT2[vars]}}{g(x)^{DT1[vars]}}_{i=1,\dots, D}
#' }
#' whereof the function g(x) stands for:
#' \deqn{g(x) = \sqrt[D]{DT[vars_1] \cdot DT[vars_2] \cdots DT[vars_D]} }
#' and "clr" is calculated by:
#' \deqn{
#'  ln \left(\frac{DT1[vars_{i}]}{DT2[vars_{i}]} * \frac{g(x)^{DT2[vars]}}{g(x)^{DT1[vars]}}\right)_{i=1,\dots, D}
#' }
#'
#' The methods "clr" and "alr" should be considered if the data contain so called \emph{compositional data} as defined by Aitchison, J. (1986): "The statistical analysis of compositional data".
#' They names correspond to the names used in the package \code{compositions} by K. Gerald van den Boogaart, Raimon Tolosana and Matevz Bren.
#'
#' Calculating the absolute error for the ratios requires calculating the absolute errors of DT1 and DT2, too.
#' For calculating the errors of DT1 and DT2 the function \code{\link{relError_dataset}} is used.
#' Accordingly the options for \code{STD_DT1} and \code{STD_DT2} are passed to the option \code{STD} in \code{relError_dataset}.
#' If STD_DT1 and/or STD_DT2 are left empty the default of 5.2\% relative error is used.
#' Also the options \code{minNr_DT1} and \code{minNr_DT2} are passed to the option \code{minNr} in \code{relError_dataset}.
#'
#' The \code{Error_method} determines how the absolute error of the ratios is calculated.
#' The error method "gauss" refers to the error propagation after Gauss:
#' \deqn{
#' \Delta x =  \frac{\Delta DT1}{DT2} - DT1 * \frac{\Delta DT2}{DT2^2}
#' }
#' The error method "biggest" refers to the maximum error after Gauss:
#' \deqn{
#' \Delta x =  \frac{\Delta DT1}{DT2} + DT1 * \frac{\Delta DT2}{DT2^2}
#' }
#'
#' For example:
#' If you have in DT1 plant samples with \code{group1.vars = "Location"} the error function would calculate the relative standard deviation for all plants of one location.
#' But maybe you have very different plants in one location so setting \code{var_subgroup = "Species"} the error function will calculate the relative standard deviation for each plant species per location, if there are more species per location than given in \code{minNr_DT1}.
#' Suppose DT2 are soil data with several samples per location.
#' If \code{group1.vars = "Location"} than the function calls \code{\link{preparation.DT2}} and calculates a mean for each location from the data set.
#' The ratio from plant to soil and the absolute errors of the ratios is then calculated for each plant sample to a mean of soils from one location.
#'
#' @return The function returns either a data.table, data.frame or a list controlled by the option \code{return_as_list}.
#' If \code{return_as_list} to FALSE a data.frame (or data.table if DT1 is of class data.table) is returned.
#' If option \code{Errors} is set to TRUE ratios and error are combined into one object and a column \code{type_of_data} is appended with the entries \emph{ratio} and \emph{ratio_error} respectively.
#' If \code{return_as_list} to TRUE the DT1-DT2-ratios are named in the list as "ratios" and, if \code{Errors} is set to TRUE the absolute errors of the ratios are saved in the list as "ratios_error".
#' If 'return_all' is set to TRUE a list with the following entries will be returned:
#'
#' [[1]] "DT1", [[2]] "DT2", [[3]] "vars", [[4]] "ratios" and if \code{Errors} is set to TRUE additionally [[5]] "ratios_error", [[6]] "DT1_error", [[7]] "DT2_error".
#'
#'
#' @family ratio functions
#' @export

ratio.DT1_DT2 <- function(DT1,
                          DT2,
                          vars = NULL,
                          group1.vars,
                          group2.vars = NULL,
                          ratio_type = "simple",
                          vars.ref,
                          Errors = FALSE,
                          Error_method = "gauss",
                          var_subgroup = NULL,
                          use_only_DT2 = FALSE,
                          DT2_replace = NULL,
                          STD_DT1,
                          STD_DT2,
                          minNr_DT1 = 50,
                          minNr_DT2 = 50,
                          return_all = FALSE,
                          return_as_list = FALSE)
{
  UseMethod("ratio.DT1_DT2", DT1)
}

#' @export
ratio.DT1_DT2.data.table <- function(DT1,
                                     DT2,
                                     vars = NULL,
                                     group1.vars,
                                     group2.vars = NULL,
                                     ratio_type = "simple",
                                     vars.ref,
                                     Errors = FALSE,
                                     Error_method = "gauss",
                                     var_subgroup = NULL,
                                     use_only_DT2 = FALSE,
                                     DT2_replace = NULL,
                                     STD_DT1,
                                     STD_DT2,
                                     minNr_DT1 = 50,
                                     minNr_DT2 = 50,
                                     return_all = FALSE,
                                     return_as_list = FALSE)
{
  if(missing(DT2)){
    print.noquote("BE AWARE: no DT2 has been provided. UpperCrust is used as DT2!")
    DT2 = data.table(UpperCrust)
  }else
    if(!is.data.table(DT2)){
      if(!is.data.frame(DT2)) stop("DT2 must be data.frame or data.table")
      DT2 = data.table(DT2)
    }

  if(missing(group1.vars)){stop("Please provide 'group1.vars'.")}
  if(!all(group1.vars %in% names(DT1))) stop(paste("In the data set of the DT1 the column(s)", group1.vars, "is missing"))
  if(!all(group1.vars %in% names(DT2))) stop(paste("In the data set of the DT2 the column(s)", group1.vars, "is missing"))
  for(kol in group1.vars){
    if(!class(DT1[[kol]]) == "factor" & !class(DT1[[kol]]) == "character") stop(paste("Column", kol, "in DT1 must be of class 'factor' or 'character'."))
    if(!class(DT2[[kol]]) == "factor" & !class(DT2[[kol]]) == "character") stop(paste("Column", kol, "in DT2 must be of class 'factor' or 'character'."))
  }
  if(!is.null(group2.vars)){
    if(length(group2.vars) > 1) stop("'group2.vars' must contain only one column name, but length is longer than 1.")
    if(!group2.vars %in% names(DT1)) stop(paste("In the data set of the DT1 the column", group2.vars, "is missing"))
    if(!class(DT1[[group2.vars]]) == "factor" & !class(DT1[[group2.vars]]) == "character") stop("Column given in 'group2.vars' in DT1 must be of class 'factor' or 'character'.")
    if(!group2.vars %in% names(DT2)) stop(paste("In the data set of the DT2 the column", group2.vars, "is missing"))
    if(!class(DT2[[group2.vars]]) == "factor" & !class(DT2[[group2.vars]]) == "character") stop("Column given in 'group2.vars' in DT2 must be of class 'factor' or 'character'.")
  }
  if(!is.null(var_subgroup)){
    if(!is.character(var_subgroup)) stop("'var_subgroup' must be a character vector.")
    if(length(var_subgroup)>1) stop("'var_subgroup' must contain only one column name, but length is longer than 1.")
    if(!var_subgroup %in% names(DT1)) stop(paste("In the data set of the DT1 the column", group2.vars, "is missing"))
    if(!class(DT1[[var_subgroup]]) == "factor" & !class(DT1[[var_subgroup]]) == "character") stop("Column given in 'var_subgroup' must be of class 'factor' or 'character'.")
  }
  # check if there are congruenting entries in column \code{group1.vars} in DT1 and DT2:
  if(sum(levels(as.factor(droplevels(DT1)[[group1.vars]])) %in% levels(as.factor(droplevels(DT2)[[group1.vars]]))) == 0){
    print.noquote(paste("In the column", group1.vars, "DT1 and DT2 have no entries in common."))
    print.noquote("Hence it is impossible to know which observations from DT2 belong to which observations in DT1.")
    if(use_only_DT2) stop("Please provide as 'group1.vars' a column which exist in both, DT1 and DT2, AND has entries in common.")
    else
      print.noquote("WARNING! DT2 is ignored and only DT2_replace is used!")
    print.noquote("")
  }

  if(is.null(vars)){
    vars = select.VarsElements(DT1, DT2)
    if(length(vars) == 0) stop("please provide for 'vars' a character vector with column names which exist in both DT1 and DT2")
  }else{
    vars_DT1 = select.VarsElements(vars, DT1)
    if(length(vars) != length(vars_DT1)){
      print.noquote("WARNING: Following 'vars' are not in DT1 and going to be ignored:")
      print.noquote(vars[!vars %in% vars_DT1])
      print.noquote("")
      vars = vars[vars %in% vars_DT1]
    }
    vars_DT2 = select.VarsElements(vars, DT2)
    if(length(vars) != length(vars_DT2)){
      print.noquote("WARNING: Following 'vars' are not in DT2 and going to be ignored:")
      print.noquote(vars[!vars %in% vars_DT2])
      print.noquote("")
      vars = vars[vars %in% vars_DT2]
    }
  }
  if(length(vars) == 0) stop("please provide for 'vars' a character vector with column names which exist in both DT1 and DT2")
  colsinfo = names(DT1)[!names(DT1) %in% vars]

  ratio_type = check_readline(ratio_type, c("simple", "ar", "alr", "cr", "clr", "log")) # check correct entry
  # make row-numbering for reordering at the end:
  DT1[, "Rows" := c(1:nrow(DT1))]


  if(Errors){ # Should errors get calculated?
    if(missing(STD_DT1)){ # Generate a data set which will make 5% error
      STD_DT1 = data.table(matrix(ncol = length(vars), nrow = 3))
      setnames(STD_DT1, vars)
      for(kol in vars){set(STD_DT1, j = kol, value = c(103.5, 100, 96.5))}
    }
    if(missing(STD_DT2)){ # Generate a data set which will make 5% error
      STD_DT2 = data.table(matrix(ncol = length(vars), nrow = 3))
      setnames(STD_DT2, vars)
      for(kol in vars){set(STD_DT2, j = kol, value = c(103.5, 100, 96.5))}
    }
    if(!all(vars %in% names(STD_DT1))){print.noquote("BE AWARE: Some errors of DT1 can't get calculated because not all columns are in STD_DT1.
                                                     Please provide a STD_DT1 which contains all 'vars', or - if you only calculate ratios and don't need the errors set 'Errors = FALSE'.")}
    if(!all(vars %in% names(STD_DT2))){print.noquote("BE AWARE: Some errors of DT2 can't get calculated because not all columns are in STD_DT2.
                                                     Please provide a STD_DT2 which contains all 'vars', or - if you only calculate ratios and don't need the errors set 'Errors = FALSE'.")}
    }else{
      STD_DT2 = NULL
      minNr_DT2 = NULL
    }

  # keying for congruenting sorting
  setkeyv(DT1, group1.vars)
  setkeyv(DT2, group1.vars)

  # prepare second data set for ratios
  DT2_prepared = preparation.DT2(DT1 = DT1, DT2 = DT2, vars = vars,
                                 group1.vars = group1.vars, group2.vars = group2.vars,
                                 use_only_DT2 = use_only_DT2, DT2_replace = DT2_replace,
                                 Errors = Errors,
                                 STD = STD_DT2, minNr = minNr_DT2, return_as_list = TRUE)

  DT2_prep = DT2_prepared$DT2 # get only the data-set
  if(use_only_DT2){ # make DT1 the same size:
    mySamples = c(rowSums(DT2_prep[, sapply(.SD, is.na), .SDcols = vars]) != length(vars)) # empty entries..
    DT1 = DT1[mySamples]
    DT2_prep = DT2_prep[mySamples]
    setkeyv(DT1, group1.vars)
    setkeyv(DT2_prep, group1.vars)
  }

  print.noquote("")
  print.noquote("-- Calculate Ratios --")
  psr = copy(DT1)
  if(ratio_type == "simple"){
    for(j in vars){set(psr, j = j, value = DT1[[j]]/DT2_prep[[j]])}
    for(j in vars){set(psr, j = j, i = which(psr[[j]]<= 0), value = NA)}
  }else{
    if(sum(is.na(DT1[, (vars), with=F]))>0){ # does DT1 contain NA in vars?
      NAColumns = colSums(DT1[, sapply(.SD, is.na), .SDcols = vars])
      if(sum(NAColumns==0) < 0.75 *length(NAColumns)){ # more than 25% of the columns have contain at least for one row NAs
        stop("There are too many vars in DT1 with NAs. Please make sure you have more columns without NAs.")
      }else{
        print.noquote("WARNING! In DT1 there are columns with NAs. If you would like to keep these columns remove the rows which contain the NAs.")
        print.noquote("Following columns are getting removed:")
        print.noquote(names(NAColumns)[NAColumns>0])
        vars = vars[!vars %in% names(NAColumns)[NAColumns>0]]
      }
    }
    if(sum(is.na(DT2_prep[, (vars), with=F]))>0){ # does DT2_prep contain NA in vars?
      NAColumns = colSums(DT2_prep[, sapply(.SD, is.na), .SDcols = vars])
      if(sum(NAColumns==0) < 0.75 *length(NAColumns)){ # more than 25% of the columns have contain at least for one row NAs
        stop("There are too many vars in DT2 with NAs. Please make sure you have more columns without NAs.")
      }else{
        print.noquote("WARNING! In DT2 there are columns with NAs. If you would like to keep these columns remove the rows which contain the NAs.")
        print.noquote("Following columns are getting removed:")
        print.noquote(names(NAColumns)[NAColumns>0])
        vars = vars[!vars %in% names(NAColumns)[NAColumns>0]]
      }
    }
  }
  if(ratio_type == "log"){
    # MeanDT1 = DT1[, apply(.SD, 1, mean, na.rm = T), .SDcols=vars]
    # MeanDT2 = DT2_prep[, apply(.SD, 1, mean, na.rm = T), .SDcols=vars]
    # meanProd = MeanDT2/MeanDT1
    for(j in vars){set(psr, j = j, value = log(DT1[[j]]/DT2_prep[[j]]))}
  }
  if(ratio_type == "cr"){
    MeanDT1 = abs(DT1[, apply(.SD, 1, prod, na.rm = T), .SDcols=vars])^c(1/length(vars))
    MeanDT2 = abs(DT2_prep[, apply(.SD, 1, prod, na.rm = T), .SDcols=vars])^c(1/length(vars))
    meanProd = MeanDT2/MeanDT1
    for(j in vars){set(psr, j = j, value = c(DT1[[j]]/DT2_prep[[j]]* meanProd))}
  }
  if(ratio_type == "clr"){
    MeanDT1 = abs(DT1[, apply(.SD, 1, prod, na.rm = T), .SDcols=vars])^c(1/length(vars))
    MeanDT2 = abs(DT2_prep[, apply(.SD, 1, prod, na.rm = T), .SDcols=vars])^c(1/length(vars))
    meanProd = MeanDT2/MeanDT1
    for(j in vars){set(psr, j = j, value = log(abs(DT1[[j]])/abs(DT2_prep[[j]])* meanProd))}
  }
  if(ratio_type == "ar"){
    if(missing(vars.ref)) vars.ref = readline(paste("Please choose one column of", paste(vars, collapse = ", "), "for calculating ar-ratios:   "))
    vars.ref = check_readline(vars.ref, vars)
    meanProd = DT2_prep[[vars.ref]]/DT1[[vars.ref]]
    for(j in vars[-which(vars %in% vars.ref)]){set(psr, j = j, value = DT1[[j]]/DT2_prep[[j]]* meanProd)}
    psr[,(vars.ref):=NULL]
  }
  if(ratio_type == "alr"){
    if(missing(vars.ref)) vars.ref = readline(paste("Please choose one column of", paste(vars, collapse = ", "), "for calculating alr-ratios:   "))
    vars.ref = check_readline(vars.ref, vars)
    meanProd = DT2_prep[[vars.ref]]/DT1[[vars.ref]]
    for(j in vars[-which(vars %in% vars.ref)]){set(psr, j = j, value = log(DT1[[j]]/DT2_prep[[j]]* meanProd))}
    psr[,(vars.ref):=NULL]
  }
  print.noquote("Done")

  if(Errors){
    print.noquote("")
    print.noquote("-- Calculate error of ratio by error of DT1 and DT2_prep --")
    DT2_abs_error = DT2_prepared$error
    if(use_only_DT2) DT2_abs_error = DT2_abs_error[mySamples]
    setkeyv(DT2_abs_error, group1.vars)

    print.noquote("")
    print.noquote("Calculate error of DT1")
    DT1_rel_error = relError_dataset(Data = DT1, vars = vars, STD = STD_DT1, group1.vars = c(var_subgroup, group1.vars),
                                     minNr = minNr_DT1)
    DT1_abs_error = copy(DT1[, vars, with = F])
    for(kol in vars){set(DT1_abs_error, j = kol, value = DT1[[kol]]*DT1_rel_error[[kol]])}

    print.noquote("")
    print.noquote("-- Calculate ratio error --")

    Error_method = check_readline(Error_method, c("gauss", "biggest"))
    if(Error_method == "gauss"){ # error through error propagation
      ratio_abs_error = abs(DT1_abs_error/DT2_prep[, vars, with = F] - (DT1[, vars, with = F]*DT2_abs_error[, vars, with = F])/DT2_prep[, vars, with = F]^2)
    }
    if(Error_method == "biggest"){ # to find biggest error the absolute value of both error-terms is used
      ratio_abs_error = DT1_abs_error/DT2_prep[, vars, with = F] + (DT1[, vars, with = F]*DT2_abs_error[, vars, with = F])/DT2_prep[, vars, with = F]^2
    }

    ratio_abs_error[ratio_abs_error<0]<-NA
    # add Infos
    ratio_abs_error = cbind(DT1[, colsinfo, with = F], abs(ratio_abs_error))
    ratio_abs_error = ratio_abs_error[rowSums(ratio_abs_error[, sapply(.SD, is.na), .SDcols = vars])!= length(vars)]
    print.noquote("Done")
  }

  if(return_all){
    if(Errors) return(list("DT1" = DT1, "DT2" = DT2_prep, "vars" = vars, "ratios" = psr, "ratios_error" = ratio_abs_error, "DT1_error" = DT1_abs_error, "DT2_error" = DT2_abs_error))
    else return(list("DT1" = DT1, "DT2" = DT2_prep, "vars" = vars, "ratios" = psr))
  }else{
    if(Errors){
      if(return_as_list) return(list("ratios" = psr, "ratios_error" = ratio_abs_error))
      else return(rbindlist(list(psr[,"type_of_data":="ratio"], ratio_abs_error[,"type_of_data":="ratio_error"]), use.names = T, fill = T))
    }else{
      if(return_as_list) return(list("ratios" = psr))
      else return(psr)
    }
  }
}

#' @export
ratio.DT1_DT2.data.frame <- function(DT1,
                                     DT2,
                                     vars = NULL,
                                     group1.vars,
                                     group2.vars = NULL,
                                     ratio_type = "simple",
                                     vars.ref,
                                     Errors = FALSE,
                                     Error_method = "gauss",
                                     var_subgroup = NULL,
                                     use_only_DT2 = FALSE,
                                     DT2_replace = NULL,
                                     STD_DT1,
                                     STD_DT2,
                                     minNr_DT1 = 50,
                                     minNr_DT2 = 50,
                                     return_all = FALSE,
                                     return_as_list = FALSE)
{
  DT1 = data.table(DT1)
  myReturn = ratio.DT1_DT2.data.table(DT1 = DT1, DT2 = DT2, vars = vars, use_only_DT2 = use_only_DT2,
                          DT2_replace = DT2_replace, group1.vars = group1.vars, group2.vars = group2.vars,
                          var_subgroup = var_subgroup, ratio_type = ratio_type, return_all = return_all,
                          return_as_list = return_as_list, Errors = Errors, Error_method = Error_method,
                          STD_DT1 = STD_DT1, STD_DT2 = STD_DT2, minNr_DT1 = minNr_DT1, minNr_DT2 = minNr_DT2,
                          vars.ref = vars.ref)
  if(return_all){
    if(Errors) return(list("DT1" = data.frame(myReturn$DT1),
                           "DT2" = data.frame(myReturn$DT2),
                           "vars" = myReturn$vars,
                           "ratios" = data.frame(myReturn$ratios),
                           "ratios_error" = data.frame(myReturn$ratios_error),
                           "DT1_error" = data.frame(myReturn$DT1_error),
                           "DT2_error" = data.frame(myReturn$DT2_error)))
    else return(list("DT1" = data.frame(myReturn$DT1),
                     "DT2" = data.frame(myReturn$DT2),
                     "vars" = myReturn$vars,
                     "ratios" = data.frame(myReturn$ratios)))
  }else{
    if(Errors){
      if(return_as_list) return(list("ratios" = data.frame(myReturn$ratios),
                                     "ratios_error" = data.frame(myReturn$ratios_error)))
      else return(data.frame(myReturn))
    }else{
      if(return_as_list) return(list("ratios" = data.frame(myReturn$ratios)))
      else return(data.frame(myReturn))
    }
  }
}

#' @export
ratio.DT1_DT2.default <- function(DT1,
                                  DT2,
                                  vars = NULL,
                                  group1.vars,
                                  group2.vars = NULL,
                                  ratio_type = "simple",
                                  vars.ref,
                                  Errors = FALSE,
                                  Error_method = "gauss",
                                  var_subgroup = NULL,
                                  use_only_DT2 = FALSE,
                                  DT2_replace = NULL,
                                  STD_DT1,
                                  STD_DT2,
                                  minNr_DT1 = 50,
                                  minNr_DT2 = 50,
                                  return_all = FALSE,
                                  return_as_list = FALSE)
{
  stop("'DT1' must be a data.frame or a data.table object.")
}

#' ratio.append_smallest
#'
#' @author Solveig Pospiech
#'
#' @description The function appends for each row the smallest ratio DT1/DT2 in the column \emph{ratio_smallest}.
#' The name of the column which contained the smallest ratio is appended in the column \emph{ratio_smallest_Elem}.
#' This function is basically a sub-function for the function \code{\link{Correction.AdheringParticles}}.
#'
#' @param Ratios list, data.frame or data.table, which is the output after using the function \code{\link{ratio.DT1_DT2}}
#' @inheritParams preparation.DT2
#'
#' @return list with [[1]] being the data set from the input with one column added containing the smallest ratio of all variables given in \code{vars}.
#' If the input was a list with one element named "ratios_error" the returned list contains a second element [[2]] "ratios_error" also with the appended columns.
#'
#' @family ratio functions
#' @export

ratio.append_smallest <- function(Ratios,
                                  vars = NULL)
{
  Errors = F
  if(!is.data.frame(Ratios)){ # if it is a list do:
    if(is.null(vars))
      if("vars" %in% names(Ratios)){vars = Ratios$vars}
    if("ratios_error" %in% names(Ratios)){
      Ratios_error = Ratios$ratios_error
      Errors = T
    }
    Ratios = Ratios$ratios
  }

  elemente = select.VarsElements(Ratios)
  if(is.null(vars)) vars = select.VarsElements(Ratios)
  else{
    # check if vars and elements are the same:
    if(!all(elemente %in% vars)){
      Ratios[, (elemente[!elemente %in% vars]) := NULL]
    }
    if(!all(vars %in% elemente)){
      print.noquote("Following elements in 'vars' are ignored because they are not present in the ratios data set:")
      print.noquote(vars[!vars %in% elemente])
      vars = vars[vars %in% elemente]
    }
  }

  Ratios = Ratios[rowSums(is.na(Ratios[, vars, with = F]))!= length(vars)] # drop rows with only NA
  Ratios[, "ratio_smallest" := apply(.SD, 1, function(x) min(x, na.rm = T)), .SDcols = vars]
  Ratios[, "ratio_smallest_Elem" := vars[apply(.SD, 1, function(x) which.min(x))], .SDcols = vars]
  if(Errors){
    for(i in 1:nrow(Ratios_error)){
      set(Ratios_error, i = i, j = "ratio_smallest", value = Ratios_error[[Ratios[["ratio_smallest_Elem"]][i]]][i])
    }
    return(list("ratios" = Ratios, "ratios_error" = Ratios_error))
  }else return(list("ratios" = Ratios))
}
