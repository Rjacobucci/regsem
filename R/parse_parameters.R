#' Takes either a vector of parameter ids or a
#' vector of named parameters and returns a vector of parameter ids
#' @param x Parameter labels
#' @param model Lavaan model
#' @return NULL if undefined input. Else vector of parameter ids
#'
parse_parameters <- function(x, model)
{
  if (is.null(x)) {return(x)} else
  if (is.numeric(x)) {return (x)}
  else if (is.character(x)) {
    labels <- parTable(model)$label
    ids <- parTable(model)$id
    matching.ids <- which(labels %in% x)
    if (length(matching.ids)!=length(x)) {
      warning("Not all requested parameters could be matched to the lavaan model parameters.")
    }
    return(ids[matching.ids])

  } else {
    warning("Unknown class type passed to parse_parameters().")
    return(NULL)
  }
}
