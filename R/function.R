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