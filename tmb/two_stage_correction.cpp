// Second-stage geostatistical correction to the dynamical model (issue #21).
//
// Gaussian model for the (working) logit bioassay mortality z_i, with its
// variance v_i fixed:
//
//   z_i = m_i + omega(s_i) + xi(s_i, t_i) + sum_k b_k[j_k(i)] + e_i,
//   e_i ~ N(0, v_i)
//
// where m_i is the dynamical model's logit prediction (an offset, never
// updated), omega is a Matern field correcting the initial conditions, xi is
// the accumulated sum since t0 of an AR(1)-in-time, Matern-in-space field of
// annual selection anomalies eta, and the b_k are iid noise effects: in the
// final model u, per pixel-year, and p, per pixel (R/two_stage_correction.R).
// The prediction target is m + omega + xi; u and p enter predictions of new
// assays only as fresh draws.
//
// Random effects:
//   w_omega  node values of omega                        (n_nodes)
//   x        node values of xi in years t0+1, ..., T     (n_nodes_xi x n_years)
//   iid      the iid effects, every term's levels stacked (n_levels)
//
// xi(., t0) = 0 is not a parameter: it enters only as the zero that the first
// time difference is taken from. x rather than eta is the random effect so
// that each observation touches the nodes of a single year; the map
// x -> eta (eta_t = x_t - x_{t-1}) is unit lower triangular, so its Jacobian
// is 1 and the prior on eta is the prior on x.
//
// With likelihood = 0 the model is Gaussian in the random effects, so given
// the hyperparameters the Laplace approximation is exact. With likelihood = 1
// (#34) the response is instead the bioassay counts, died_i of n_i,
// beta-binomial with mean ilogit(lambda_i) and intra-class correlation rho_i,
// and TMB's Laplace approximation integrates the random effects out of that
// likelihood itself; z and v are then unused.

#include <TMB.hpp>

using namespace density;
using namespace R_inla;
using namespace Eigen;

// log prior density of a d = 2 Matern field's (log kappa, log sigma), the
// scale the optimiser works on, with range r = sqrt(8) / kappa: log r ~
// normal(prior(0), prior(1)) (|d log r / d log kappa| = 1), and the Fuglstad
// et al. (2019) PC prior on sigma, P(sigma > prior(2)) = prior(3)
// (|d sigma / d log sigma| = sigma)
template<class Type>
Type log_matern_prior(Type log_kappa, Type log_sigma, vector<Type> prior) {
  Type log_range = log(sqrt(Type(8.0))) - log_kappa;
  Type sigma = exp(log_sigma);
  Type lambda_sigma = -log(prior(3)) / prior(2);
  Type lp = dnorm(log_range, prior(0), prior(1), true);
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
  DATA_VECTOR(z);                 // (working) logit response
  DATA_VECTOR(v);                 // its fixed variance
  DATA_VECTOR(m);                 // dynamical model logit prediction (offset)
  DATA_SPARSE_MATRIX(A_omega);    // observations -> omega mesh nodes
  DATA_SPARSE_MATRIX(A_xi);       // observations -> vec(x) (node x year)
  DATA_IMATRIX(iid_index);        // 0-based level in iid, observation x term
  DATA_IVECTOR(iid_term);         // 0-based term of each level in iid
  DATA_STRUCT(spde, spde_t);      // FEM matrices of the omega mesh
  DATA_STRUCT(spde_xi, spde_t);   // and of the xi mesh
  DATA_VECTOR(prior_omega);       // (mean log range, sd log range, sigma0, P(sigma > sigma0))
  DATA_VECTOR(prior_eta);         // as above, for the eta innovations field
  DATA_VECTOR(pc_iid);            // (sigma0, P(sigma > sigma0)) for every iid term
  DATA_VECTOR(persistence_prior); // (meanlog, sdlog) of 1 / (1 - phi)
  DATA_INTEGER(likelihood);       // 0: Gaussian working response, 1: beta-binomial counts
  DATA_VECTOR(died);              // the counts, for likelihood = 1
  DATA_VECTOR(n_mosquito);
  DATA_VECTOR(rho);               // beta-binomial intra-class correlation

  PARAMETER_VECTOR(w_omega);
  PARAMETER_ARRAY(x);
  PARAMETER_VECTOR(iid);
  PARAMETER(log_sigma_omega);
  PARAMETER(log_kappa_omega);
  PARAMETER(log_sigma_eta);
  PARAMETER(log_kappa_eta);
  PARAMETER(logit_phi);
  PARAMETER_VECTOR(log_sigma_iid); // one per iid term

  Type phi = invlogit(logit_phi);
  vector<Type> sigma_iid = exp(log_sigma_iid);
  Type nll = 0.0;

  // omega: Matern field at the nodes
  nll += GMRF(matern_precision(spde, exp(log_kappa_omega),
                               exp(log_sigma_omega)))(w_omega);
  vector<Type> lambda = m + A_omega * w_omega;

  // xi: accumulated AR(1) of Matern fields. eta holds the time differences of
  // x (with x_t0 = 0). SEPARABLE(f, g) applies f to the last (year) dimension
  // and g to the first (node); TMB's AR1 is the unit-variance stationary
  // process, so eta has stationary marginal SD sigma_eta, from the GMRF
  int n_nodes = x.dim(0);
  int n_years = x.dim(1);
  array<Type> eta(n_nodes, n_years);
  for (int i = 0; i < n_nodes; i++) {
    eta(i, 0) = x(i, 0);
    for (int t = 1; t < n_years; t++) eta(i, t) = x(i, t) - x(i, t - 1);
  }
  nll += SEPARABLE(AR1(phi), GMRF(matern_precision(spde_xi,
                                                   exp(log_kappa_eta),
                                                   exp(log_sigma_eta))))(eta);
  vector<Type> x_vec = x.vec();
  lambda += A_xi * x_vec;

  // iid effects
  for (int j = 0; j < iid.size(); j++) {
    nll -= dnorm(iid(j), Type(0.0), sigma_iid(iid_term(j)), true);
  }
  for (int i = 0; i < z.size(); i++) {
    for (int k = 0; k < iid_index.cols(); k++) lambda(i) += iid(iid_index(i, k));
  }

  if (likelihood == 0) {
    nll -= dnorm(z, lambda, sqrt(v), true).sum();
  } else {
    // log p and log (1 - p) without overflow; the binomial coefficient is a
    // constant and left out
    for (int i = 0; i < z.size(); i++) {
      Type log_p = -logspace_add(Type(0.0), -lambda(i));
      Type log_q = -logspace_add(Type(0.0), lambda(i));
      Type c = (Type(1.0) - rho(i)) / rho(i);
      Type a = c * exp(log_p);
      Type b = c * exp(log_q);
      nll -= lgamma(died(i) + a) + lgamma(n_mosquito(i) - died(i) + b) -
        lgamma(n_mosquito(i) + c) - lgamma(a) - lgamma(b) + lgamma(c);
    }
  }

  // priors on the hyperparameters: PC priors on both Matern fields,
  // exponential (PC) priors on the iid SDs (with the Jacobian for log sigma),
  // and a lognormal prior on persistence L = 1 / (1 - phi) = 1 + exp(logit_phi),
  // for which d log L / d logit_phi = phi
  nll -= log_matern_prior(log_kappa_omega, log_sigma_omega, prior_omega);
  nll -= log_matern_prior(log_kappa_eta, log_sigma_eta, prior_eta);
  Type lambda_iid = -log(pc_iid(1)) / pc_iid(0);
  for (int k = 0; k < sigma_iid.size(); k++) {
    nll -= log(lambda_iid) - lambda_iid * sigma_iid(k) + log_sigma_iid(k);
  }
  Type persistence = Type(1.0) / (Type(1.0) - phi);
  nll -= dnorm(log(persistence), persistence_prior(0), persistence_prior(1),
               true) + log(phi);

  return nll;
}
