#' Plot Method for a `summary_fitsae` Object
#'
#' The generic method `plot()` provides, in a grid (default) or sequence, (a) a scatterplot of direct estimates versus model-based estimates, visually capturing the shrinking process, (b) a Bayesian P-values histogram, (c) a boxplot of standard deviation reduction values, and, if areas sample sizes are provided as input in `fit_sae()`, (d) a scatterplot of model residuals versus sample sizes, in order to check for design-consistency i.e., as long as sizes increase residuals should converge to zero.
#'
#' @param x Object of class `summary_fitsae`.
#' @param size Aesthetic option denoting the size of scatterplots points, see \code{\link[ggplot2]{geom_point}} documentation.
#' @param alpha Aesthetic option denoting the opacity of scatterplots points, see \code{\link[ggplot2]{geom_point}} documentation.
#' @param n_bins Denoting the number of bins used for histogram.
#' @param grid Logical indicating whether plots are displayed in a grid (`TRUE`) or in sequence (`FALSE`).
#' @param label_names Character string indicating the model name to display in boxplot x-axis label.
#' @param ... Currently unused.
#'
#' @return Four `ggplot2` objects in a grid.
#'
#' @seealso \code{\link{summary.fitsae}} to produce the input object.
#' @examples
#' library(tipsae)
#'
#' # loading toy dataset
#' data("emilia_cs")
#'
#' # fitting a model
#' fit_beta <- fit_sae(formula_fixed = hcr ~ x, data = emilia_cs, domains = "id",
#'                     type_disp = "var", disp_direct = "vars", domain_size = "n",
#'                     # MCMC setting to obtain a fast example. Remove next line for reliable results.
#'                     chains = 1, iter = 150, seed = 0)
#'
#' # check model diagnostics
#' summ_beta <- summary(fit_beta)
#'
#' # visualize diagnostics via plot() method
#' plot(summ_beta)
#'
#' @export

plot.summary_fitsae <- function(x,
                   size = 2.5,
                   alpha = 0.8,
                   n_bins = 15,
                   grid = TRUE,
                   label_names = NULL,
                   ...
                   ){
  if (!inherits(x, "summary_fitsae"))
    stop("Indicated object does not have 'summary_fitsae' class.")

  if (is.null(x$data_obj$domain_size_n))
    message("When checking for design consistency, the domains sample size is required. \n
            Specify the argument 'domain_size' of the fit_sae() function.")

  if (!grid) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))}

    # Arranging dataset
    xydata <- data.frame(x = x$direct_est,
                         y = x$post_means,
                         res = x$residuals,
                         bp = x$bayes_pvalues,
                         ind = ifelse(!is.null(label_names), label_names, "Mod1")
                         )
    if (is.null(x$data_obj$domain_size_n)) {
      xydata$n <- NA
    }else{
      xydata$n <- x$data_obj$domain_size_n
    }
    if (is.null(x$sd_reduction)) {
      xydata$sdr <- NA
    }else{
      xydata$sdr <- x$sd_reduction
    }

    # Boxplot for standard deviation reduction
    boxplot_sdr <- ggplot2::ggplot(xydata, ggplot2::aes_(y = ~ sdr, x = ~ ind)) +
      ggplot2::theme_bw() + ggplot2::ylab("S.D. Reduction") + ggplot2::xlab("Distribution") +
      ggplot2::geom_boxplot() +
      ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank())
    if (0 >= min(xydata$sdr) & 0 <= max(xydata$sdr))
      boxplot_sdr <- boxplot_sdr + ggplot2::geom_hline(yintercept = 0)

    # Plot direct vs model estimates
    lims_axis <- range(c(x$direct_est, x$post_means))
    scatter_s <- ggplot2::ggplot(data = xydata, ggplot2::aes_(x = ~ x, y = ~ y)) +
      ggplot2::geom_abline(slope = 1, intercept = 0) +
      ggplot2::xlim(lims_axis) + ggplot2::ylim(lims_axis) +
      ggplot2::theme(aspect.ratio = 1) +
      ggplot2::ylab("HB est.") +
      ggplot2::xlab("Direct est.") +
      ggplot2::theme_bw() +
      ggplot2::geom_point(
        shape = 20,
        size = size,
        alpha = alpha
      )

    # Histogram Bayesian p-values
    hist_bp <- ggplot2::ggplot(data = xydata, ggplot2::aes_(x = ~ bp)) +
      ggplot2::xlim(0, 1) + ggplot2::theme_bw() +
      ggplot2::xlab("Bayesian p-values") + ggplot2::ylab("Counts") +
      ggplot2::geom_histogram(bins = n_bins,
                              color = "black",
                              fill = "white")

    # Design consistency
    if (!is.null(x$data_obj$domain_size_n)) {
      scatter_dc <- ggplot2::ggplot(data = xydata, ggplot2::aes_(x = ~ n, y = ~ res)) +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::theme_bw() +
        ggplot2::geom_point(
          shape = 20,
          size = size,
          alpha = alpha
        ) + ggplot2::ylab("Residuals") +
        ggplot2::xlab("Domain sample size")

    if (grid) {
      gridExtra::grid.arrange(scatter_s, hist_bp,
                              boxplot_sdr, scatter_dc,
                              ncol = 2)
    }else{
    print(scatter_s)
    print(hist_bp)
    print(boxplot_sdr)
    print(scatter_dc)
    }
  }else{
    if (grid) {
      gridExtra::grid.arrange(scatter_s, hist_bp,
                              boxplot_sdr,
                              ncol = 2)
      }else{
        print(scatter_s)
        print(hist_bp)
        print(boxplot_sdr)
      }
  }
}
