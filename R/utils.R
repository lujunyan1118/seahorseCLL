# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


#' Create pretty scientific notation
#'
#' Function to generate pretty scientific notation format for plot label
#' @param n Number for formating
#' @param digits Digits of decimals to keep
#' @export
sciPretty <- function(n, digits = 2, bold = FALSE) {
  nForm <- strsplit(format(n, digits = digits, scientific = TRUE),split = "e")
  b <- nForm[[1]][1]
  i <- as.integer(nForm[[1]][2])
  #bquote(.(b)%*%10^.(i))
  if(bold) {
    sprintf("bold(%s%%*%%10^%s)",b,i)
  } else sprintf("%s%%*%%10^%s",b,i)
}

#' Format Seahorse names
#'
#' Function to give nice names of seahorse measurement
#' @param seaList A character of Seahorse measurement names
#' @export

formatSea <- function(seaList) {
  newList <- sapply(seaList, function(x) gsub("[.]"," ",x))
  newList[newList == "ECAR OCR ratio"] <- "ECAR/OCR"
  return(newList)
}

#' Calculate modified row z-scores
#'
#' Row centered by median and scaled by median absolute deviation
#' @param x The input matrix or datafram
#' @param center Centered by median
#' @param scale Scaled by 1.4826*MAD
#' @export

mscale <- function(x, center = TRUE, scale = TRUE){
  x.scaled <- apply(x, 1, function(y) (y-median(y, na.rm = TRUE))/(1.4826*mad(y, na.rm = TRUE)))
  return(t(as.matrix(x.scaled)))
}

