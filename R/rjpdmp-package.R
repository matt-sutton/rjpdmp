#'rjpdmp-package
#'
#' @description
#' Implements various reversible jump piecewise-deterministic Markov Process methods including the ZigZag and Bouncy Particle Sampler with Normal or Spherical velocity distributions (Chevallier, Fearnhead, Sutton 2020, https://arxiv.org/abs/2010.11771).
#'
#' The package can be used to Generates PDMP trajectories for reversible jump
#' \itemize{
#'      \item \code{\link{zigzag_logit}}: ZigZag on logistic likelihood problem
#'      \item \code{\link{zigzag_logit_ss}}: ZigZag with subsampling on Logistic likelihood problem
#'      \item \code{\link{bps_s_logit}}: BPS with velocities distributed uniformly on the sphere for a Logistic likelihood problem
#'      \item \code{\link{bps_n_logit}}: BPS with velocities distributed Normally for a Logistic likelihood problem
#'      \item \code{\link{zigzag_rr}}: ZigZag on a robust regression likelihood problem
#'      \item \code{\link{bps_s_rr}}: BPS with velocities distributed uniformly on the sphere for a robust regression likelihood problem
#'      \item \code{\link{bps_n_rr}}: BPS with velocities distributed Normally for a robust regression likelihood problem
#'}
#'
#'@section Additional functions:
#'
#' Additional functions for plotting, generating samples, calculating posterior means or probabilities of inclusion
#' \itemize{
  #'      \item \code{\link{plot_pdmp}}: Plot marginal densities and joint pairs plots for trajectories and samples of PDMP samplers and optionally MCMC samples for comparison.
#'      \item \code{\link{plot_pdmp_multiple}}: Plots to compare PDMP samplers and optionally MCMC samples.
  #'      \item \code{\link{gen_sample}}: Get samples from PDMP trajectories taking a fixed time discretisation.
#'      \item \code{\link{model_probabilities}}: Calculate either marginal probabilities of inclusions or posterior probabilities of specific models.
  #'      \item \code{\link{models_visited}}: Count the number of times a model is visited
#'      \item \code{\link{marginal_mean}}: Calculate the marginal mean using PDMP trajectories
  #'      \item \code{\link{cond_mean}}: Calculate the mean conditioned on being in a specific model
#'      }
#'
#' Extensions to the package are planned.
#'
"_PACKAGE"
#> [1] "_PACKAGE"
