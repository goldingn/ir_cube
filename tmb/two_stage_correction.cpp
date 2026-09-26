// Second-stage geostatistical correction to the dynamical model (issue #21).
//
// Gaussian model for the empirical-logit bioassay mortality z_i, with its
// sampling variance v_i fixed:
//
//   z_i = m_i + omega(s_i) + xi(s_i, t_i) + u_j(i) [+ p_c(i)] [+ s_k(i)] + e_i,
//   e_i ~ N(0, v_i)
//
// where m_i is the dynamical model's logit prediction (an offset, never
// updated), omega is a Matern field correcting the initial conditions, xi is
// the accumulated sum since t0 of an AR(1)-in-time, Matern-in-space field of
// annual selection anomalies eta, and u is an iid pixel-year effect.
//
// Two optional iid terms extend the error structure (doc/two_stage_plan.md,
// "Error structure"); both are off by default, which gives exactly the model
// without them:
//   p  static pixel effect, one per pixel (raster cell), SD sigma_p: a
//      persistent local deviation, part of the prediction target;
//   s  survey effect, one per survey (citation x country x year), SD sigma_s:
//      shared measurement / batch error, not part of the prediction target.
// When a term is off its vector is mapped to zero in R, and neither its prior
// nor its hyperprior enters the objective.
//
// Random effects:
//   w_omega  node values of omega                           (n_nodes)
//   x        node values of xi in years t0+1, ..., T        (n_nodes_xi x n_years)
//
// xi may use a coarser mesh than omega: its latent dimension is multiplied by
// the number of years, and the fill-in of the Cholesky factor of the
// space-time block dominates the cost of the fit.
//   u        pixel-year effects                             (n_pixel_years)
//   p        static pixel effects                           (n_pixels, or 1 if off)
//   s        survey effects                                 (n_surveys, or 1 if off)
// xi(., t0) = 0 is not a parameter: it enters only as the zero that the first
// time difference is taken from.
//
// x rather than eta is the random effect so that each observation touches the
// nodes of a single year: parameterising by eta would make every observation
// depend on all earlier years' eta, and the Hessian much denser. The map
// x -> eta (first differences, with x_t0 = 0) is unit lower triangular, so its
// Jacobian is 1 and the prior on eta is the prior on x.
//
// The model is Gaussian in the random effects, so given the hyperparameters the
// Laplace approximation is exact and the inner Newton step converges at once.

#include <TMB.hpp>

using namespace density;
using namespace R_inla;
using namespace Eigen;

// log density of the Fuglstad et al. (2019) PC prior for a d = 2 Matern field,
// on the (log kappa, log sigma) scale the optimiser works on.
//
// With range r = sqrt(8) / kappa, the prior on (r, sigma) is
//   pi(r)     = lambda_r r^-2 exp(-lambda_r / r),  lambda_r = -log(alpha_r) r0
//   pi(sigma) = lambda_s exp(-lambda_s sigma),    lambda_s = -log(alpha_s) / s0
// so that P(r < r0) = alpha_r and P(sigma > s0) = alpha_s. The Jacobians are
// |d r / d log kappa| = r and |d sigma / d log sigma| = sigma.
template<class Type>
Type log_pc_matern(Type log_kappa, Type log_sigma, vector<Type> pc) {
  Type range0 = pc(0);
  Type alpha_range = pc(1);
  Type sigma0 = pc(2);
  Type alpha_sigma = pc(3);
  Type log_range = log(sqrt(Type(8.0))) - log_kappa;
  Type range = exp(log_range);
  Type sigma = exp(log_sigma);
  Type lambda_range = -log(alpha_range) * range0;
  Type lambda_sigma = -log(alpha_sigma) / sigma0;
  Type lp = log(lambda_range) - Type(2.0) * log_range - lambda_range / range +
    log_range;
  lp += log(lambda_sigma) - lambda_sigma * sigma + log_sigma;
  return lp;
}

// SPDE precision (alpha = 2, d = 2) scaled to have marginal SD sigma
template<class Type>
SparseMatrix<Type> matern_precision(spde_t<Type> spde, Type kappa,
                                    Type sigma) {
  Type tau_spde = Type(1.0) / (sigma * kappa * sqrt(Type(4.0 * M_PI)));
  SparseMatrix<Type> Q = Q_spde(spde, kappa);
  return Q * (tau_spde * tau_spde);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // data -----------------------------------------------------------------
  DATA_VECTOR(z);                 // empirical logit
  DATA_VECTOR(v);                 // its (fixed) sampling variance
  DATA_VECTOR(m);                 // dynamical model logit prediction (offset)
  DATA_SPARSE_MATRIX(A_omega);    // observations -> mesh nodes
  DATA_SPARSE_MATRIX(A_xi);       // observations -> vec(x) (node x year)
  DATA_IVECTOR(u_index);          // 0-based pixel-year of each observation
  DATA_STRUCT(spde, spde_t);      // FEM matrices M0 (c0), M1 (g1), M2 (g2)
  DATA_STRUCT(spde_xi, spde_t);   // as above, for the (possibly coarser) xi mesh
  DATA_INTEGER(include_xi);       // 1 for variant omega_xi_u, 0 for omega_u
  DATA_VECTOR(pc_omega);          // (range0, P(range < range0), sigma0, P(sigma > sigma0))
  DATA_VECTOR(pc_eta);            // as above, for the eta innovations field
  DATA_VECTOR(pc_tau);            // (tau0, P(tau > tau0))
  DATA_VECTOR(persistence_prior); // (meanlog, sdlog) of 1 / (1 - phi)
  DATA_INTEGER(include_p);        // 1 to include the static pixel effect p
  DATA_IVECTOR(p_index);          // 0-based pixel of each observation
  DATA_VECTOR(pc_sigma_p);        // (sigma0, P(sigma_p > sigma0))
  DATA_INTEGER(include_s);        // 1 to include the survey effect s
  DATA_IVECTOR(s_index);          // 0-based survey of each observation
  DATA_VECTOR(pc_sigma_s);        // (sigma0, P(sigma_s > sigma0))

  // parameters -------------------------------------------------------------
  PARAMETER_VECTOR(w_omega);
  PARAMETER_ARRAY(x);
  PARAMETER_VECTOR(u);
  PARAMETER(log_sigma_omega);
  PARAMETER(log_kappa_omega);
  PARAMETER(log_sigma_eta);
  PARAMETER(log_kappa_eta);
  PARAMETER(logit_phi);
  PARAMETER(log_tau);
  PARAMETER_VECTOR(p);
  PARAMETER_VECTOR(s);
  PARAMETER(log_sigma_p);
  PARAMETER(log_sigma_s);

  Type sigma_omega = exp(log_sigma_omega);
  Type kappa_omega = exp(log_kappa_omega);
  Type sigma_eta = exp(log_sigma_eta);
  Type kappa_eta = exp(log_kappa_eta);
  Type phi = invlogit(logit_phi);
  Type tau = exp(log_tau);

  Type nll = 0.0;

  // omega: Matern field at the nodes ----------------------------------------
  SparseMatrix<Type> Q_omega = matern_precision(spde, kappa_omega, sigma_omega);
  nll += GMRF(Q_omega)(w_omega);

  vector<Type> lambda = m + A_omega * w_omega;

  // xi: integrated AR(1) of Matern fields -------------------------------------
  if (include_xi == 1) {
    int n_nodes = x.dim(0);
    int n_years = x.dim(1);

    // time differences, with x_t0 = 0
    array<Type> eta(n_nodes, n_years);
    for (int i = 0; i < n_nodes; i++) {
      eta(i, 0) = x(i, 0);
      for (int t = 1; t < n_years; t++) {
        eta(i, t) = x(i, t) - x(i, t - 1);
      }
    }

    // SEPARABLE(f, g) applies f to the last (outer, slowest-varying) array
    // dimension and g to the first, so with eta stored nodes x years this is
    // AR(1) over years and GMRF over nodes. TMB's AR1 is the unit-variance
    // stationary process eta_t = phi eta_{t-1} + sqrt(1 - phi^2) eps_t with
    // eta_1 ~ N(0, 1): exactly the issue's scaling, so eta has stationary
    // marginal SD sigma_eta, taken from the GMRF.
    SparseMatrix<Type> Q_eta = matern_precision(spde_xi, kappa_eta, sigma_eta);
    nll += SEPARABLE(AR1(phi), GMRF(Q_eta))(eta);

    vector<Type> x_vec = x.vec();
    lambda += A_xi * x_vec;
  }

  // u: pixel-year effects -----------------------------------------------------
  nll -= dnorm(u, Type(0.0), tau, true).sum();
  for (int i = 0; i < z.size(); i++) {
    lambda(i) += u(u_index(i));
  }

  // p: static pixel effects, s: survey effects (optional) -------------------
  Type sigma_p = exp(log_sigma_p);
  Type sigma_s = exp(log_sigma_s);
  if (include_p == 1) {
    nll -= dnorm(p, Type(0.0), sigma_p, true).sum();
    for (int i = 0; i < z.size(); i++) {
      lambda(i) += p(p_index(i));
    }
  }
  if (include_s == 1) {
    nll -= dnorm(s, Type(0.0), sigma_s, true).sum();
    for (int i = 0; i < z.size(); i++) {
      lambda(i) += s(s_index(i));
    }
  }

  // likelihood ------------------------------------------------------------------
  nll -= dnorm(z, lambda, sqrt(v), true).sum();

  // priors (penalties) on the hyperparameters ----------------------------------
  nll -= log_pc_matern(log_kappa_omega, log_sigma_omega, pc_omega);

  // exponential (PC) prior on tau, with the Jacobian for log tau
  Type lambda_tau = -log(pc_tau(1)) / pc_tau(0);
  nll -= log(lambda_tau) - lambda_tau * tau + log_tau;

  // exponential (PC) priors on sigma_p and sigma_s, as for tau
  if (include_p == 1) {
    Type lambda_p = -log(pc_sigma_p(1)) / pc_sigma_p(0);
    nll -= log(lambda_p) - lambda_p * sigma_p + log_sigma_p;
  }
  if (include_s == 1) {
    Type lambda_s = -log(pc_sigma_s(1)) / pc_sigma_s(0);
    nll -= log(lambda_s) - lambda_s * sigma_s + log_sigma_s;
  }

  Type persistence = Type(1.0) / (Type(1.0) - phi);
  if (include_xi == 1) {
    nll -= log_pc_matern(log_kappa_eta, log_sigma_eta, pc_eta);
    // lognormal prior on persistence L = 1 / (1 - phi) = 1 + exp(logit_phi),
    // so d log L / d logit_phi = phi
    nll -= dnorm(log(persistence), persistence_prior(0), persistence_prior(1),
                 true) + log(phi);
  }

  // reports ---------------------------------------------------------------------
  Type range_omega = sqrt(Type(8.0)) / kappa_omega;
  Type range_eta = sqrt(Type(8.0)) / kappa_eta;
  REPORT(range_omega);
  REPORT(sigma_omega);
  REPORT(range_eta);
  REPORT(sigma_eta);
  REPORT(phi);
  REPORT(persistence);
  REPORT(tau);
  REPORT(sigma_p);
  REPORT(sigma_s);
  ADREPORT(range_omega);
  ADREPORT(sigma_omega);
  ADREPORT(tau);
  if (include_xi == 1) {
    ADREPORT(range_eta);
    ADREPORT(sigma_eta);
    ADREPORT(phi);
  }

  return nll;
}
