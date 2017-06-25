#' Correction.AdheringParticles
#'
#' @author Solveig Pospiech
#'
#' @description Suppose element data of one data set (DT1) are biased because the concentrations are the result of a mixture of two substances, of which one substance are the element concentrations of DT2.
#' In order to correct DT1 to \eqn{DT_{corrected}} a fraction of DT2 has to be subtracted from DT1.
#' The basic equation for the correction is:
#' \deqn{
#' DT_{corrected}=\frac{DT1 - x * DT2}{1 - x}
#' }
#' whereof \emph{x} is the amount of DT2 to be subtracted.
#'
#' The function is written for the case that \emph{x} is unknown.
#' To calculate \emph{x} the condition is that in \eqn{DT_{corrected}} at least one element concentration is zero or known.
#' Suppose \eqn{vars_{i}} has a very low concentration, close to zero, in \eqn{DT_{corrected}}: \eqn{DT_{corrected}[vars_{i}]=0}, then:
#' \deqn{
#' x =  \frac{DT1[vars_{i}]}{DT2[vars_{i}]}
#' }
#'
#' The function was developed for the use to correct plant concentrations for adhering particles:
#' Exact and reproducible analysis of element concentrations in plant tissue is the basis for many research fields such as environmental, health, phytomining, agricultural or provenance studies.
#' Unfortunately plant samples collected in the field will always contain particles on their tissue surfaces such as airborne dust or soil particles.
#' If not removed these particles may induce a bias to the element concentrations measured in plant samples.
#'
#' For full description of the calculations and the background of correction plants for adhering particles please refer to:
#'
#' Pospiech, S., Fahlbusch, W., Sauer, B., Pasold, T., & Ruppert, H. (2017). Alteration of trace element concentrations in plants by adhering particles–Methods of correction. \emph{Chemosphere, 182, 501-508}.
#' and the section Details.
#'
#' @inheritParams preparation.DT2
#' @inheritParams ratio.DT1_DT2
#' @param vars_ignore character vector of column names, only for 'method 3'.
#' These variables are ignored for calculating the median of amount of DT2 (\emph{x}) in 'method 3'.
#' Please note: the functions returns corrected values for these columns because they are only ignored for calculating the median of \emph{x}.
#' Default is "As", "Se", "Sn", "V", "Be", "Ge" and "Pt". Please see Details for further explanation.
#' @param id.vars column with unique (!) entries for each row. Class can be integer (corresponding row numbers) or character (e.g. sample IDs).
#' If missing, all columns but \code{vars} will be assigned to it.
#' Please note: Function is faster and more stable if \code{id.vars} is provided.
#' @param method characters (no character vector!, please give m3 instead of "m3") denoting the method. Options are \emph{m1}, \emph{m2} and \emph{m3} and \emph{subtr}.
#' Default is m3. Please see details.
#' @param element string, only for method 1. Denotes the column with which amount of DT2 (\emph{x}) is to be calculated.
#' @param offset numeric, default is 0. The offset diminishes the subtracted amount of DT2 \emph{x}: x = x - offset.
#' If used with m2 all concentrations will stay > 0.
#' Reasonable offset is e.g. offset = 0.0001
#' @param negative_values logical, should negative values be returned? If set to FALSE negative values are set to 0. Default is FALSE.
#' @param set_statistical_0 logical, only for method 3.
#' Should all values of the variables contributing to the median of \emph{x} be set to 0? Default is FALSE.
#'
#' @details The main option of this function is the \code{method} which determines how the amount of DT2 to be subtracted, the \emph{x}, is going to be calculated.
#' There are four options:
#' \itemize{
#' \item Method 1: calculate \emph{x} via a fixed element
#' \item Method 2: calculate \emph{x} via the element with the smallest ratio between DT1[vars] and DT2[vars]
#' \item Method 3: calculate \emph{x} via the median of several, very small ratios between DT1[vars] and DT2[vars]
#' \item Method subtr: calculate the concentrations for \eqn{x * DT2[vars]}
#' }
#'
#' To Method 1:
#' For example using Ti as \code{element} \eqn{DT_{corrected}} is calculated with \eqn{ x = DT1[Ti]/DT2[Ti]}.
#' Typical elements for the option \code{element} are e.g. Ti, Al, Zr, Sc, ...
#' This will eventually lead to negative concentrations for some elements.
#'
#' To Method 2:
#' This method subtracts the smallest possible content of DT2 from DT1 (smallest \emph{x}).
#' For each row/sample the element with the smallest \emph{x} of all ratios \eqn{ x = DT1[vars]/DT2[vars]} of each sample is taken as \emph{element},
#' hence every sample is corrected based on a different \emph{element}.
#' With this method there are no negative concentrations.
#'
#' To Method 3:
#' In order to reduce the uncertainty of the content of DT2 in DT1 (\emph{x}) based on only one element as in method 1 and 2 an average of several \emph{x} of elements can be calculated.
#' With \eqn{\Delta x} being the absolute error of \emph{x} the median is calculated by all \emph{x} of elements which values \eqn{ x - \Delta x} are smaller than \eqn{ x_{smallest} + \Delta x_{smallest}}.
#' The value of the median \eqn{\bar{x}} is then used as \emph{x}.
#' This will eventually lead to negative concentrations for some elements.
#' Because statistically the \emph{x} of all elements, which error overlaps the error of the element with smallest \emph{x}, are indistinguishable we suggest to set all elements contributing to \eqn{\bar{x}} to zero, because these small values should not be interpreted:
#' Set option \code{set_statistical_0} to TRUE.
#'
#' It is advisable to exclude elements with a huge error margin in the option \code{vars_ignore} because they could severely increase the median \eqn{\bar{x}} by "opening" the window of error-ranges for many elements with significantly higher ratios.
#' This could lead to an unnatural high median \eqn{\bar{x}} resulting into an overcorrection.
#'
#' If option \code{id.vars} is provided the functions prints the 'group1.vars' and 'id.vars' of the sample.
#'
#' For examples and more information please refer to:
#' Pospiech, S., Fahlbusch, W., Sauer, B., Pasold, T., & Ruppert, H. (2017). Alteration of trace element concentrations in plants by adhering particles–Methods of correction. \emph{Chemosphere, 182, 501-508}.
#'
#'
#' @return data.frame (or data.table if DT1 is data.table) according to \code{method}.
#'
#' @family ratio functions
#' @export

Correction.AdheringParticles<- function(DT1,
                                        DT2 = NULL,
                                        vars = NULL,
                                        vars_ignore = c("As", "Se", "Sn", "V", "Be", "Ge", "Pt"),
                                        method,
                                        element,
                                        id.vars,
                                        group1.vars,
                                        group2.vars,
                                        var_subgroup,
                                        offset = 0,
                                        use_only_DT2 = TRUE,
                                        DT2_replace = NULL,
                                        Errors = TRUE,
                                        return_as_list = TRUE,
                                        negative_values = FALSE,
                                        set_statistical_0 = FALSE,
                                        Error_method = "gauss",
                                        STD_DT1 = STD_Plant,
                                        STD_DT2 = STD_Soil,
                                        minNr_DT1 = 100,
                                        minNr_DT2 = 100)
{
  UseMethod("Correction.AdheringParticles", DT1)
}

#' @export
Correction.AdheringParticles.data.table <- function(DT1,
                                                    DT2 = NULL,
                                                    vars = NULL,
                                                    vars_ignore = c("As", "Se", "Sn", "V", "Be", "Ge", "Pt"),
                                                    method,
                                                    element,
                                                    id.vars,
                                                    group1.vars,
                                                    group2.vars,
                                                    var_subgroup,
                                                    offset = 0,
                                                    use_only_DT2 = TRUE,
                                                    DT2_replace = NULL,
                                                    Errors = TRUE,
                                                    return_as_list = TRUE,
                                                    negative_values = FALSE,
                                                    set_statistical_0 = FALSE,
                                                    Error_method = "gauss",
                                                    STD_DT1 = STD_Plant,
                                                    STD_DT2 = STD_Soil,
                                                    minNr_DT1 = 100,
                                                    minNr_DT2 = 100)
{
  if(missing(group1.vars)){stop("Please provide 'group1.vars'.")}
  if(!is.character(group1.vars)) stop("'group1.vars' must be a character vector.")
  if(!group1.vars %in% names(DT1)) stop(paste("In the data set of the DT1 the column", group1.vars, "is missing"))
  if(missing(DT2)){
    print.noquote("BE AWARE: no DT2 has been provided. UpperCrust is used as DT2!")
    DT2 = data.table(UpperCrust)
  }else
    if(!is.data.table(DT2)){
      if(!is.data.frame(DT2)) stop("DT2 must be data.frame or data.table")
      DT2 = data.table(DT2)
    }
  if(missing(group2.vars)){group2.vars = NULL}
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
  if(missing(var_subgroup)){var_subgroup = NULL}

  if(missing(method)){
    print.noquote("Please choose a method for calculation:")
    print.noquote("There are two methods available - M1 (M2) and M3.")
    print.noquote("All are described in the Paper 'Alteration of trace element concentrations in
                  DT1 by adhering particles - methods of correction'.")
    print.noquote("M1 and M2 are in this function the same because which Element is used is determined in the function 'ratio.append_smallest'.")
    print.noquote("M3 is method 3 as described in the paper.")
    print.noquote("'dust' calculates the amount of dust on the leaves.")
    method = readline("Which method you would like to use?")
  }
  method = tolower(method)
  method = str_trim(method)
  method = check_readline(method, c("m1", "m2", "m3", "subtracted"))

  if(method == "m3"){ # calculate necessary with errors:
    Ratios_list = ratio.DT1_DT2(DT1 = DT1, DT2 = DT2, vars = vars, use_only_DT2 = use_only_DT2, DT2_replace = DT2_replace,
                                group1.vars = group1.vars, group2.vars = group2.vars, STD_DT1 = STD_DT1, STD_DT2 = STD_DT2,
                                return_all = T, minNr_DT1 = minNr_DT1, minNr_DT2 = minNr_DT2, Errors = TRUE,
                                var_subgroup = var_subgroup, Error_method = Error_method, return_as_list = TRUE)
  }else{ # leave Errors to option
    Ratios_list = ratio.DT1_DT2(DT1 = DT1, DT2 = DT2, vars = vars, use_only_DT2 = use_only_DT2, DT2_replace = DT2_replace,
                                group1.vars = group1.vars, group2.vars = group2.vars, STD_DT1 = STD_DT1, STD_DT2 = STD_DT2,
                                return_all = T, minNr_DT1 = minNr_DT1, minNr_DT2 = minNr_DT2, Errors = Errors,
                                var_subgroup = var_subgroup, Error_method = Error_method, return_as_list = TRUE)
  }

  DT1 = Ratios_list$DT1
  vars = Ratios_list$vars
  cols_taken = vars[!vars %in% vars_ignore]

  if(method == c("m1")){
    ratios = Ratios_list$ratios
    ratios_error = Ratios_list$ratios_error
    if(missing(element)){
      element = readline("Correction with which fixed element?   ")}
    ratios[, "ratio_smallest" := get(element)]
    ratios[, "ratio_smallest_Elem"  := element]
    ratios_error[, "ratio_smallest" := get(element)]
    ratios_error[, "ratio_smallest_Elem"  := element]
  }else{
    Ratios = ratio.append_smallest(Ratios = Ratios_list, vars = cols_taken)
    ratios = Ratios$ratios
    ratios_error = Ratios$ratios_error
  }

  DT1_abs_error = Ratios_list$DT1_error
  DT1_abs_error[, (group1.vars) := DT1[[group1.vars]]]
  # DT2_prep
  DT2_prep = Ratios_list$DT2
  DT2_prep_abs_error = Ratios_list$DT2_error

  # add offset
  set(ratios, j = "ratio_smallest", value = ratios$ratio_smallest - offset)
  set(ratios_error, j = "ratio_smallest", i = which(ratios[["ratio_smallest"]] < 0), value = 0)
  set(ratios, j = "ratio_smallest", i = which(ratios[["ratio_smallest"]] < 0), value = 0)

  setkeyv(ratios, group1.vars)
  setkeyv(ratios_error, group1.vars)
  setkeyv(DT1, group1.vars)
  setkeyv(DT2_prep, group1.vars)
  setkeyv(DT2_prep_abs_error, group1.vars)
  setkeyv(DT1_abs_error, group1.vars)

  # because for values below detection limit the values might be marked as "-xx" all negative values are set to zero:
  for(j in vars){set(DT1, j = j, i = which(DT1[[j]]<0), value = NA)}
  for(j in vars){set(DT2_prep, j = j, i = which(DT2_prep[[j]]<0), value = NA)}
  for(j in cols_taken){set(ratios, j = j, i = which(ratios[[j]]<0), value = NA)}

  # ----- Dustcorrection ----

  # m2 ----
  if(method == "m2"){method = "m1"}

  # m1 ----
  if(method == "m1"){
    sk_DT1 = (DT1[, vars, with = F]-(ratios$ratio_smallest*DT2_prep[, vars, with = F]))/(1-ratios$ratio_smallest)
    sk_DT1[, paste(vars) := round(.SD, digits = 10), .SDcols = vars]

    suppressWarnings(sk_DT1[sk_DT1<0]<-0) # set all negtive values to zero

    # absolute error for DT1
    DT1_abs_error = sqrt(
      (DT1_abs_error[, vars, with = F]/(1-ratios$ratio_smallest))^2 +
        ((ratios$ratio_smallest*DT2_prep_abs_error[, vars, with = F])/(1-ratios$ratio_smallest))^2  +
        (matrix(ratios_error$ratio_smallest/(1-ratios$ratio_smallest)^2, ncol = length(vars), nrow = nrow(sk_DT1)))^2
    )
  }

  if(method == "m3"){ # m3 ----
    # Bei welchen Proben und Elementen überlappen sich die Fehler?
    TrueFalseMatrix = data.table(matrix(nrow = nrow(ratios), ncol = length(cols_taken)))
    setnames(TrueFalseMatrix, cols_taken)
    for(j in cols_taken){
      set(TrueFalseMatrix, j = j,
          value = c(!is.na(ratios[[j]]) & (ratios[[j]] - ratios_error[[j]]) <= (ratios[["ratio_smallest"]] + ratios_error[["ratio_smallest"]])))
    }
    TrueFalseMatrix[is.na(TrueFalseMatrix)]<-FALSE

    print.noquote("")
    print.noquote("Which elements are how often contributing to the ratio final:")
    print.noquote(sort(colSums(TrueFalseMatrix)))
    print.noquote("")

    # neue ratio_smallest bei den TF und bei den Fehlern für output error:
    zw = ratios[, cols_taken, with = F]
    zw[!TrueFalseMatrix]<-NA
    if(sum(rowSums(zw[, sapply(.SD, is.na), .SDcols = cols_taken]) == length(cols_taken))>0){
      print.noquote("following samples are not changed due to offset:")
      if(missing(id.vars)){
        print(DT1[rowSums(zw[, sapply(.SD, is.na), .SDcols = cols_taken]) == length(cols_taken), group1.vars, with = F])
      }else{
        print(DT1[rowSums(zw[, sapply(.SD, is.na), .SDcols = cols_taken]) == length(cols_taken), c(id.vars, group1.vars), with = F])
      }
    }
    ratio_Final = apply(zw, 1, function(x) median(x, na.rm = T)) # median aus allen SF, die im Bereich des Fehlers liegen
    ratio_Final[is.na(ratio_Final)]<-0

    zw = ratios_error[, cols_taken, with = F]
    zw[!TrueFalseMatrix]<-NA
    ratio_FinalFehler = apply(zw, 1, function(x) median(x, na.rm = T)) # median aus allen SF, die im Bereich des Fehlers liegen
    ratio_FinalFehler = ratio_FinalFehler[is.na(ratio_FinalFehler)]<-0

    sk_DT1 = (DT1[, vars, with = F]-(ratio_Final*DT2_prep[, vars, with = F]))/(1-ratio_Final)
    sk_DT1[, paste(vars) := round(.SD, digits = 10), .SDcols = vars]

    if(set_statistical_0){ # values to 0 if there are statistically not distinguishable from 0
      for(j in cols_taken){ set(sk_DT1, j = j, i = which(TrueFalseMatrix[[j]]), value = 0)}
    }

    # absolute error for corrected DT1
    if(Errors){
      DT1_abs_error = sqrt(
        (DT1_abs_error[, vars, with = F]/(1-ratio_Final))^2 +
          ((ratio_Final*DT2_prep_abs_error[, vars, with = F])/(1-ratio_Final))^2  +
          (matrix(ratio_FinalFehler/(1-ratio_Final)^2, ncol = length(vars), nrow = nrow(sk_DT1)))^2
      )
      ratios_error$ratio_smallest = ratio_FinalFehler
    }
    ratios$ratio_smallest = ratio_Final
  }

  if(method == "subtracted"){
    DT1[, (vars) := ratios$ratio_smallest*DT2_prep/(1-ratios$ratio_smallest)]
    myData = DT1
    return(myData)
  }else{
    if(!negative_values){ # set negative values to 0, if wished
      sk_DT1[sk_DT1<0]<-0
    }
    for(kol in names(DT1)[!names(DT1) %in% vars]){set(sk_DT1, j = kol, value = DT1[[kol]])} # add all other columns than vars to sk_DT1
    sk_DT1[, "psr" := ratios$ratio_smallest] # add the 'psr'-column for knowing the calculated amount of adhering particles
    if(Errors){
      DT1_corr_error = data.table(DT1_abs_error)
      for(kol in names(DT1)[!names(DT1) %in% vars]){set(DT1_corr_error, j = kol, value = DT1[[kol]])} # add all other columns than vars to errors
      DT1_corr_error[, "psr" := ratios_error$ratio_smallest]
    }

    #Return data
    if(Errors){
      if(return_as_list){
        sk_DT1 = sk_DT1[rowSums(sk_DT1[, sapply(.SD, is.na), .SDcols = vars]) != length(vars)]
        setkeyv(sk_DT1, group1.vars)
        DT1_corr_error = DT1_corr_error[rowSums(DT1_corr_error[, sapply(.SD, is.na), .SDcols = vars]) != length(vars)]
        setkeyv(DT1_corr_error, group1.vars)
        return(list(DT1_corr = sk_DT1, Errors = DT1_corr_error))
      }else{
        sk_DT1[, "KindOfResult" := "conc"]
        DT1_corr_error[, "KindOfResult" := "error"]
        myData = rbindlist(list(sk_DT1, DT1_corr_error)) # combine as data.table
        myData = myData[rowSums(myData[, sapply(.SD, is.na), .SDcols = vars]) != length(vars)]
        setkeyv(myData, group1.vars)
        return(myData)
      }
    }else{
      myData = sk_DT1[rowSums(sk_DT1[, sapply(.SD, is.na), .SDcols = vars]) != length(vars)]
      setkeyv(myData, group1.vars)
      return(myData)
    }
  }
}

#' @export
Correction.AdheringParticles.data.frame <- function(DT1,
                                                    DT2 = NULL,
                                                    vars = NULL,
                                                    vars_ignore = c("As", "Se", "Sn", "V", "Be", "Ge", "Pt"),
                                                    method,
                                                    element,
                                                    id.vars,
                                                    group1.vars,
                                                    group2.vars,
                                                    var_subgroup,
                                                    offset = 0,
                                                    use_only_DT2 = TRUE,
                                                    DT2_replace = NULL,
                                                    Errors = TRUE,
                                                    return_as_list = TRUE,
                                                    negative_values = FALSE,
                                                    set_statistical_0 = FALSE,
                                                    Error_method = "gauss",
                                                    STD_DT1 = STD_Plant,
                                                    STD_DT2 = STD_Soil,
                                                    minNr_DT1 = 100,
                                                    minNr_DT2 = 100)
{
  DT1 = data.table(DT1)
  myReturn = Correction.AdheringParticles.data.table(DT1 = DT1, DT2 = DT2, vars = vars, vars_ignore = vars_ignore,
                                                     method = method, id.vars = id.vars, group1.vars = group1.vars,
                                                     group2.vars = group2.vars, var_subgroup = var_subgroup, offset = offset,
                                                     use_only_DT2 = use_only_DT2, DT2_replace = DT2_replace, Errors = Errors,
                                                     return_as_list = return_as_list, negative_values = negative_values,
                                                     set_statistical_0 = set_statistical_0, Error_method = Error_method,
                                                     minNr_DT1 = minNr_DT1, minNr_DT2 = minNr_DT2)
  #Return data
  if(Errors){
    if(return_as_list){
      return(list("DT1_corr" = data.frame(myReturn$DT1_corr), "Errors" = data.frame(myReturn$Errors)))
    }else{
      return(data.frame(myReturn))
    }
  }else{
    return(data.frame(myReturn))
  }
}

#' @export
Correction.AdheringParticles.default <- function(DT1,
                                                 DT2 = NULL,
                                                 vars = NULL,
                                                 vars_ignore = c("As", "Se", "Sn", "V", "Be", "Ge", "Pt"),
                                                 method,
                                                 element,
                                                 id.vars,
                                                 group1.vars,
                                                 group2.vars,
                                                 var_subgroup,
                                                 offset = 0,
                                                 use_only_DT2 = TRUE,
                                                 DT2_replace = NULL,
                                                 Errors = TRUE,
                                                 return_as_list = TRUE,
                                                 negative_values = FALSE,
                                                 set_statistical_0 = FALSE,
                                                 Error_method = "gauss",
                                                 STD_DT1 = STD_Plant,
                                                 STD_DT2 = STD_Soil,
                                                 minNr_DT1 = 100,
                                                 minNr_DT2 = 100)
{
  stop("'DT1' must be data.frame or data.table")
}
