# functions

orchard_plot2 <- function (object, mod = "Int", xlab, N = "none", alpha = 0.5, 
          angle = 90, cb = FALSE, k = TRUE, transfm = c("none", "tanh"), point.size = 3,
          condition.lab = "Condition") 
{
  transfm <- match.arg(transfm)
  if (any(class(object) %in% c("rma.mv", "rma"))) {
    if (mod != "Int") {
      object <- mod_results(object, mod)
    }
    else {
      object <- mod_results(object, mod = "Int")
    }
  }
  mod_table <- object$mod_table
  data <- object$data
  data$moderator <- factor(data$moderator, levels = mod_table$name, 
                           labels = mod_table$name)
  data$scale <- (1/sqrt(data[, "vi"]))
  legend <- "Precision (1/SE)"
  if (any(N != "none")) {
    data$scale <- N
    legend <- "Sample Size (N)"
  }
  if (transfm == "tanh") {
    cols <- sapply(mod_table, is.numeric)
    mod_table[, cols] <- Zr_to_r(mod_table[, cols])
    data$yi <- Zr_to_r(data$yi)
    label <- xlab
  }
  else {
    label <- xlab
  }
  mod_table$K <- as.vector(by(data, data[, "moderator"], function(x) length(x[, 
                                                                              "yi"])))
  group_no <- length(unique(mod_table[, "name"]))
  cbpl <- c("#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
            "#CC79A7", "#56B4E9", "#999999")
  if (names(mod_table)[2] == "condition") {
    condition_no <- length(unique(mod_table[, "condition"]))
    plot <- ggplot2::ggplot() + ggbeeswarm::geom_quasirandom(data = data, 
                                                             ggplot2::aes(y = yi, x = moderator, size = scale, 
                                                                          colour = moderator), alpha = alpha) + ggplot2::geom_hline(yintercept = 0, 
                                                                                                                                    linetype = 2, colour = "black", alpha = alpha) + 
      ggplot2::geom_linerange(data = mod_table, ggplot2::aes(x = name, 
                                                             ymin = lowerCL, ymax = upperCL), size = 1.2, 
                              position = ggplot2::position_dodge2(width = 0.3)) + 
      ggplot2::geom_pointrange(data = mod_table, ggplot2::aes(y = estimate, 
                                                              x = name, ymin = lowerPR, ymax = upperPR, shape = as.factor(condition), 
                                                              fill = name), size = 0.5, position = ggplot2::position_dodge2(width = 0.3)) + 
      ggplot2::scale_shape_manual(values = 20 + (1:condition_no)) + 
      ggplot2::coord_flip() + ggplot2::theme_bw() + ggplot2::guides(fill = "none", 
                                                                    colour = "none") + ggplot2::theme(legend.position = c(0, 
                                                                                                                          1), legend.justification = c(0, 1)) + ggplot2::theme(legend.title = ggplot2::element_text(size = 9)) + 
      ggplot2::theme(legend.direction = "horizontal") + 
      ggplot2::theme(legend.background = ggplot2::element_blank()) + 
      ggplot2::labs(y = label, x = "", size = legend) + 
      ggplot2::labs(shape = condition.lab) + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, 
                                                                                                colour = "black", hjust = 0.5, angle = angle))
    plot <- plot + ggplot2::annotate("text", y = (max(data$yi) + 
                                                    (max(data$yi) * 0.1)), x = (seq(1, group_no, 1) + 
                                                                                  0.3), label = paste("italic(k)==", mod_table$K[1:group_no]), 
                                     parse = TRUE, hjust = "right", size = 3.5)
  }
  else {
    plot <- ggplot2::ggplot(data = mod_table, ggplot2::aes(x = estimate, 
                                                           y = name)) + ggbeeswarm::geom_quasirandom(data = data, 
                                                                                                     ggplot2::aes(x = yi, y = moderator, size = scale, 
                                                                                                                  colour = moderator), groupOnX = FALSE, alpha = alpha) + 
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = lowerPR, 
                                           xmax = upperPR), height = 0, show.legend = FALSE, 
                              size = 0.5, alpha = 0.6) + ggplot2::geom_errorbarh(ggplot2::aes(xmin = lowerCL, 
                                                                                              xmax = upperCL), height = 0, show.legend = FALSE, 
                                                                                 size = 1.2) + ggplot2::geom_vline(xintercept = 0, 
                                                                                                                   linetype = 2, colour = "black", alpha = alpha) + 
      ggplot2::geom_point(ggplot2::aes(fill = name), size = point.size, 
                          shape = 21) + ggplot2::theme_bw() + ggplot2::guides(fill = "none", 
                                                                              colour = "none") + ggplot2::theme(legend.position = c(1, 
                                                                                                                                    0), legend.justification = c(1, 0)) + ggplot2::theme(legend.title = ggplot2::element_text(size = 9)) + 
      ggplot2::theme(legend.direction = "horizontal") + 
      ggplot2::theme(legend.background = ggplot2::element_blank()) + 
      ggplot2::labs(x = label, y = "", size = legend) + 
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10, 
                                                         colour = "black", hjust = 0.5, angle = angle))
    if (k == TRUE) {
      plot <- plot + ggplot2::annotate("text", x = (max(data$yi) + 
                                                      (max(data$yi) * 0.1)), y = (seq(1, group_no, 
                                                                                      1) + 0.3), label = paste("italic(k)==", mod_table$K), 
                                       parse = TRUE, hjust = "right", size = 3.5)
    }
  }
  if (cb == TRUE) {
    plot <- plot + ggplot2::scale_fill_manual(values = cbpl) + 
      ggplot2::scale_colour_manual(values = cbpl)
  }
  return(plot)
}

# mod_result old

#' @title get_est
#' @description Function gets estimates from rma objects (metafor)
#' @param model rma.mv object
#' @param mod the name of a moderator. If meta-analysis (i.e. no moderator, se mod = "Int")
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @export

get_est <- function (model, mod) {
  name <- firstup(as.character(stringr::str_replace(row.names(model$beta), {{mod}}, "")))
  
  estimate <- as.numeric(model$beta)
  lowerCL <- model$ci.lb
  upperCL <- model$ci.ub
  
  table <- tibble::tibble(name = factor(name, levels = name, labels = name), estimate = estimate, lowerCL = lowerCL, upperCL = upperCL)
  
  return(table)
}


#' @title get_pred
#' @description Function to get prediction intervals (crediblity intervals) from rma objects (metafor)
#' @param model rma.mv object
#' @param mod the name of a moderator
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @export

get_pred <- function (model, mod) {
  
  name <- firstup(as.character(stringr::str_replace(row.names(model$beta), {{mod}}, "")))
  len <- length(name)
  
  if(len != 1){
    newdata <- matrix(NA, ncol = len, nrow = len)
    
    pred <- metafor::predict.rma(model, newmods = diag(len),
                                 tau2.levels = 1:len,
                                 gamma2.levels = 1:len)
  }
  else {
    pred <- metafor::predict.rma(model)
  }
  lowerPR <- pred$cr.lb
  upperPR <- pred$cr.ub
  
  table <- tibble::tibble(name = factor(name, levels = name, labels = name), lowerPR = lowerPR, upperPR = upperPR)
  return(table)
}

#' @title firstup
#' @description Uppercase moderator names
#' @param x a character string
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a character string with all combinations of the moderator level names with upper case first letters
#' @export
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

#' @title get_data
#' @description Collects and builds the data used to fit the rma.mv or rma model in metafor
#' @param model rma.mv object
#' @param mod the moderator variable
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a data frame
#' @export
#'
get_data <- function(model, mod){
  X <- as.data.frame(model$X)
  names <- vapply(stringr::str_split(colnames(X), {{mod}}), function(x) paste(unique(x), collapse = ""), character(1L))
  
  moderator <- matrix(ncol = 1, nrow = dim(X)[1])
  
  for(i in 1:ncol(X)){
    moderator <- ifelse(X[,i] == 1, names[i], moderator)
  }
  moderator <- firstup(moderator)
  yi <- model$yi
  vi <- model$vi
  type <- attr(model$yi, "measure")
  
  data <- data.frame(yi, vi, moderator, type)
  return(data)
  
}

#' @title mod_results
#' @description Using a metafor model object of class rma or rma.mv it creates a table of model results containing the mean effect size estimates for all levels of a given categorical moderator, their corresponding confidence intervals and prediction intervals
#' @param model rma.mv object
#' @param mod the name of a moderator; put "Int" if the intercept model (meta-analysis) or no moderators.
#' @return A data frame containing all the model results including mean effect size estimate, confidence and prediction intervals
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{data(eklof)
#' eklof<-metafor::escalc(measure="ROM", n1i=N_control, sd1i=SD_control,
#' m1i=mean_control, n2i=N_treatment, sd2i=SD_treatment, m2i=mean_treatment,
#' data=eklof)
#' # Add the unit level predictor
#' eklof$Datapoint<-as.factor(seq(1, dim(eklof)[1], 1))
#' # fit a MLMR - accouting for some non-independence
#' eklof_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~ Grazer.type-1, random=list(~1|ExptID,
#' ~1|Datapoint), data=eklof)
#' results <- mod_results(eklof_MR, mod = "Grazer.type")
#' }
#' @export

mod_results <- function(model, mod) {
  
  if(all(class(model) %in% c("rma.mv", "rma.uni", "rma")) == FALSE) {stop("Sorry, you need to fit a metafor model of class rma.mv or rma")}
  
  data <- get_data(model, mod)
  
  # Get confidence intervals
  CI <- get_est(model, mod)
  
  # Get prediction intervals
  PI <- get_pred(model, mod)
  
  model_results <- list(mod_table = cbind(CI, PI[,-1]), data = data)
  
  class(model_results) <- "orchard"
  
  return(model_results)
  
}
# TODO - I think we can improve `mod` bit?

#' @title print.orchard
#' @description Print method for class 'orchard'
#' @param object x an R object of class orchard
#' @param ... Other arguments passed to print
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @return Returns a data frame
#' @export
#'
print.orchard <- function(object, ...){
  return(object$mod_table)
}
