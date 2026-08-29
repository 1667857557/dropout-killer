#include <RcppArmadillo.h>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

static bool dk_inverse(const arma::mat& A, arma::mat& inv) {
  bool ok = arma::inv_sympd(inv, A);
  if (!ok) ok = arma::inv(inv, A);
  return ok && inv.is_finite();
}

// [[Rcpp::export]]
Rcpp::List dk_batch_ridge_positive_cpp(
    const arma::mat& y,
    const arma::mat& scores,
    const Rcpp::IntegerVector& query_gene,
    const Rcpp::IntegerVector& query_cell,
    const double ridge,
    const int min_target_observed,
    const bool support_adaptive_rank,
    const double bias_kappa) {
  const int G = y.n_rows, n = y.n_cols, E = query_gene.size();
  NumericVector prediction(E, 0.0), factor_prediction(E, NA_REAL),
    prediction_sd(E, NA_REAL), predictability(E, 0.0), shrinkage(E, 0.0),
    bias_calibration(E, 0.0);
  IntegerVector n_observed(E, 0), factor_rank_used(E, 0), method(E, 0);
  std::vector<std::vector<int> > event_by_gene(G);
  for (int e = 0; e < E; ++e) {
    int g = query_gene[e] - 1;
    if (g >= 0 && g < G) event_by_gene[g].push_back(e);
  }
  for (int g = 0; g < G; ++g) {
    const std::vector<int>& qe = event_by_gene[g];
    if (qe.empty()) continue;
    std::vector<arma::uword> donor_idx;
    donor_idx.reserve(n);
    double sum1 = 0.0, sum2 = 0.0;
    for (int j = 0; j < n; ++j) {
      const double v = y(g, j);
      if (std::isfinite(v) && v > 0.0) {
        donor_idx.push_back(j); sum1 += v; sum2 += v * v;
      }
    }
    const int nobs = donor_idx.size();
    const double mu = nobs ? sum1 / nobs : 0.0;
    const double total_sse = nobs ? std::max(0.0, sum2 - sum1 * sum1 / nobs) : 0.0;
    const double total_var = nobs > 1 ? total_sse / (nobs - 1.0) : NA_REAL;
    auto fallback = [&]() {
      for (int e : qe) {
        n_observed[e] = nobs;
        if (nobs >= 2 && std::isfinite(total_var)) {
          prediction[e] = std::max(0.0, mu);
          prediction_sd[e] = std::sqrt(std::max(0.0, total_var * (1.0 + 1.0 / nobs)));
          method[e] = 1;
        }
      }
    };
    if (nobs < 2 || !std::isfinite(total_var)) { fallback(); continue; }
    int k = scores.n_cols;
    if (support_adaptive_rank) k = std::min(k, std::max(0, (nobs - 3) / 2));
    if (k < 1 || nobs < std::max(min_target_observed, k + 3) || total_sse <= 1e-12) {
      fallback(); continue;
    }
    const int p = k + 1;
    arma::mat X(nobs, p, arma::fill::ones);
    arma::vec target(nobs);
    for (int r = 0; r < nobs; ++r) {
      const int j = donor_idx[r]; target[r] = y(g, j);
      for (int a = 0; a < k; ++a) X(r, a + 1) = scores(j, a);
    }
    arma::mat A = X.t() * X;
    for (int a = 1; a < p; ++a) A(a, a) += ridge;
    arma::mat inv;
    if (!dk_inverse(A, inv)) { fallback(); continue; }
    arma::vec beta = inv * X.t() * target;
    arma::vec fitted = X * beta;
    arma::vec hdiag = arma::sum((X * inv) % X, 1);
    arma::vec loo_factor(nobs), loo_null(nobs), d(nobs), centered(nobs);
    for (int r = 0; r < nobs; ++r) {
      double denom = std::max(1e-6, 1.0 - hdiag[r]);
      loo_factor[r] = target[r] - (target[r] - fitted[r]) / denom;
      loo_null[r] = (sum1 - target[r]) / (nobs - 1.0);
      d[r] = loo_factor[r] - loo_null[r];
      centered[r] = target[r] - loo_null[r];
    }
    double den = arma::dot(d, d);
    double q = den > 1e-12 ? arma::dot(d, centered) / den : 0.0;
    q = std::max(0.0, std::min(1.0, q));
    arma::vec loo = loo_null + q * d;
    double loo_bias = arma::mean(loo - target);
    double bias = std::isfinite(bias_kappa) ? nobs / (nobs + bias_kappa) * loo_bias : 0.0;
    if (!std::isfinite(bias)) bias = 0.0;
    arma::vec loo_cal = loo - bias;
    double null_sse = arma::dot(centered, centered);
    double model_sse = arma::dot(target - loo_cal, target - loo_cal);
    double predability = null_sse > 0 ? 1.0 - model_sse / null_sse : 0.0;
    predability = std::max(0.0, std::min(1.0, predability));
    double sigma2 = model_sse / nobs;
    if (!std::isfinite(sigma2) || sigma2 < 0) sigma2 = total_var;
    for (int e : qe) {
      int c = query_cell[e] - 1;
      arma::rowvec xq(p, arma::fill::ones);
      for (int a = 0; a < k; ++a) xq[a + 1] = scores(c, a);
      double raw = arma::as_scalar(xq * beta);
      double pred = std::max(0.0, mu + q * (raw - mu) - bias);
      double hq = arma::as_scalar(xq * inv * xq.t());
      double leverage = (1.0 - q) * (1.0 - q) / nobs + q * q * std::max(0.0, hq);
      double pv = sigma2 * (1.0 + leverage);
      if (!std::isfinite(pv) || pv < 0) pv = total_var * (1.0 + 1.0 / nobs);
      prediction[e] = pred; factor_prediction[e] = raw;
      prediction_sd[e] = std::sqrt(std::max(0.0, pv));
      predictability[e] = predability; shrinkage[e] = q;
      bias_calibration[e] = bias; n_observed[e] = nobs;
      factor_rank_used[e] = k; method[e] = q > 0 ? 2 : 1;
    }
  }
  return List::create(
    _["prediction"] = prediction,
    _["factor_prediction"] = factor_prediction,
    _["prediction_sd"] = prediction_sd,
    _["predictability"] = predictability,
    _["shrinkage"] = shrinkage,
    _["bias_calibration"] = bias_calibration,
    _["n_observed"] = n_observed,
    _["factor_rank_used"] = factor_rank_used,
    _["method"] = method
  );
}

// [[Rcpp::export]]
Rcpp::List dk_batch_tree_ridge_positive_cpp(
    const arma::mat& y,
    const arma::mat& local_mean,
    const arma::mat& local_variance,
    const arma::mat& effective_n,
    const arma::mat& W,
    const arma::mat& scores,
    const Rcpp::IntegerVector& query_gene,
    const Rcpp::IntegerVector& query_cell,
    const double ridge,
    const int min_target_observed,
    const double min_effective_donors,
    const double local_info_kappa,
    const bool support_adaptive_rank) {
  const int G = y.n_rows, n = y.n_cols, E = query_gene.size();
  NumericVector prediction(E, NA_REAL), residual_prediction(E, NA_REAL),
    prediction_sd(E, NA_REAL), predictability(E, 0.0), shrinkage(E, 0.0);
  IntegerVector n_donors(E, 0), factor_rank_used(E, 0), method(E, 0);
  std::vector<std::vector<int> > event_by_gene(G);
  for (int e = 0; e < E; ++e) {
    int g = query_gene[e] - 1;
    if (g >= 0 && g < G) event_by_gene[g].push_back(e);
  }
  for (int g = 0; g < G; ++g) {
    const std::vector<int>& qe = event_by_gene[g];
    if (qe.empty()) continue;
    std::vector<int> supported_events;
    for (int e : qe) {
      int c = query_cell[e] - 1;
      double mu = local_mean(g, c), vv = local_variance(g, c), ne = effective_n(g, c);
      if (std::isfinite(mu) && std::isfinite(ne) && ne >= min_effective_donors) {
        prediction[e] = std::max(0.0, mu); method[e] = 1;
        if (std::isfinite(vv) && vv >= 0) prediction_sd[e] = std::sqrt(vv + vv / ne);
        supported_events.push_back(e);
      }
    }
    if (supported_events.empty() || scores.n_cols < 1) continue;
    std::vector<int> donors;
    for (int j = 0; j < n; ++j)
      if (std::isfinite(y(g, j)) && y(g, j) > 0 && std::isfinite(local_mean(g, j))) donors.push_back(j);
    if (donors.empty()) continue;
    arma::vec agg(donors.size(), arma::fill::zeros);
    for (int e : supported_events) {
      int c = query_cell[e] - 1;
      for (unsigned int r = 0; r < donors.size(); ++r) agg[r] += W(c, donors[r]);
    }
    std::vector<int> kept;
    for (unsigned int r = 0; r < donors.size(); ++r)
      if (std::isfinite(agg[r]) && agg[r] > 0) kept.push_back(r);
    const int nd = kept.size();
    int k = scores.n_cols;
    if (support_adaptive_rank) k = std::min(k, std::max(0, (nd - 3) / 2));
    if (k < 1 || nd < std::max(min_target_observed, k + 2)) continue;
    const int p = k + 1;
    arma::mat X(nd, p, arma::fill::ones);
    arma::vec target(nd), weight(nd);
    for (int r = 0; r < nd; ++r) {
      int old = kept[r], j = donors[old];
      target[r] = y(g, j) - local_mean(g, j); weight[r] = agg[old];
      for (int a = 0; a < k; ++a) X(r, a + 1) = scores(j, a);
    }
    double mw = arma::mean(weight);
    if (!std::isfinite(mw) || mw <= 0) continue;
    arma::vec scale_w = weight / mw, sw = arma::sqrt(scale_w);
    arma::mat Xw = X.each_col() % sw;
    arma::vec yw = target % sw;
    arma::mat A = Xw.t() * Xw;
    for (int a = 1; a < p; ++a) A(a, a) += ridge;
    arma::mat inv;
    if (!dk_inverse(A, inv)) continue;
    arma::vec beta = inv * Xw.t() * yw;
    arma::vec fitted = X * beta;
    arma::vec hdiag = scale_w % arma::sum((X * inv) % X, 1);
    arma::vec loo(nd);
    for (int r = 0; r < nd; ++r) {
      double denom = std::max(1e-6, 1.0 - hdiag[r]);
      loo[r] = target[r] - (target[r] - fitted[r]) / denom;
    }
    double den = arma::dot(weight, loo % loo);
    double qpred = den > 1e-12 ? arma::dot(weight, loo % target) / den : 0.0;
    qpred = std::max(0.0, std::min(1.0, qpred));
    arma::mat meat = X.t() * (X.each_col() % arma::square(scale_w));
    arma::mat M = inv * meat * inv;
    for (int e : qe) {
      int c = query_cell[e] - 1;
      double mu = local_mean(g, c), vv = local_variance(g, c), ne = effective_n(g, c);
      bool supported = std::isfinite(mu) && std::isfinite(ne) && ne >= min_effective_donors;
      if (!supported) continue;
      arma::rowvec xq(p, arma::fill::ones);
      for (int a = 0; a < k; ++a) xq[a + 1] = scores(c, a);
      double raw = arma::as_scalar(xq * beta);
      double qinfo = local_info_kappa > 0 ? ne / (ne + local_info_kappa) : 1.0;
      double qfinal = qpred * qinfo;
      double pred = std::max(0.0, mu + qfinal * raw);
      double null_sse = 0, cross = 0, loo_sse = 0, wsum = 0;
      int ndonor_q = 0;
      for (int r = 0; r < nd; ++r) {
        int j = donors[kept[r]]; double w = W(c, j);
        if (w > 0) ndonor_q++;
        null_sse += w * target[r] * target[r];
        cross += w * target[r] * loo[r];
        loo_sse += w * loo[r] * loo[r]; wsum += w;
      }
      double final_sse = std::max(0.0, null_sse - 2 * qfinal * cross + qfinal * qfinal * loo_sse);
      double pa = null_sse > 0 ? 1.0 - final_sse / null_sse : 0.0;
      pa = std::max(0.0, std::min(1.0, pa));
      double sigma2 = wsum > 0 ? final_sse / wsum : vv;
      if (!std::isfinite(sigma2) || sigma2 < 0) sigma2 = vv;
      double base_var = (std::isfinite(vv) && ne > 0) ? vv / ne : 0.0;
      double hq = arma::as_scalar(xq * M * xq.t());
      double pv = sigma2 * (1.0 + qfinal * qfinal * std::max(0.0, hq)) + base_var;
      prediction[e] = pred; residual_prediction[e] = raw;
      prediction_sd[e] = std::sqrt(std::max(0.0, pv));
      predictability[e] = pa; shrinkage[e] = qfinal;
      n_donors[e] = ndonor_q; factor_rank_used[e] = k;
      method[e] = qfinal > 0 ? 2 : 1;
    }
  }
  return List::create(
    _["prediction"] = prediction,
    _["residual_prediction"] = residual_prediction,
    _["prediction_sd"] = prediction_sd,
    _["predictability"] = predictability,
    _["shrinkage"] = shrinkage,
    _["n_donors"] = n_donors,
    _["factor_rank_used"] = factor_rank_used,
    _["method"] = method
  );
}

// Build target-safe P1 predictions for every cell and analytic LOO predictions
// for positive donors.  The R architecture dispatcher consumes the state once
// for several recovery routes, avoiding repeated ridge fits.
// [[Rcpp::export]]
Rcpp::List dk_batch_p1_state_cpp(
    const arma::mat& y,
    const arma::mat& scores,
    const double ridge,
    const int min_target_observed,
    const bool support_adaptive_rank,
    const double bias_kappa) {
  const int G = y.n_rows, n = y.n_cols;
  arma::mat prediction(G, n, arma::fill::zeros);
  arma::mat prediction_sd(G, n); prediction_sd.fill(NA_REAL);
  arma::mat loo(G, n); loo.fill(NA_REAL);
  NumericVector shrinkage(G, 0.0), bias_calibration(G, 0.0), sigma2_out(G, NA_REAL);
  IntegerVector n_observed(G, 0), factor_rank_used(G, 0), method(G, 0);
  for (int g = 0; g < G; ++g) {
    std::vector<arma::uword> donor_idx;
    donor_idx.reserve(n);
    double sum1 = 0.0, sum2 = 0.0;
    for (int j = 0; j < n; ++j) {
      const double v = y(g, j);
      if (std::isfinite(v) && v > 0.0) {
        donor_idx.push_back(j); sum1 += v; sum2 += v * v;
      }
    }
    const int nobs = donor_idx.size(); n_observed[g] = nobs;
    const double mu = nobs ? sum1 / nobs : 0.0;
    const double sse = nobs ? std::max(0.0, sum2 - sum1 * sum1 / nobs) : 0.0;
    const double total_var = nobs > 1 ? sse / (nobs - 1.0) : NA_REAL;
    auto fallback = [&]() {
      prediction.row(g).fill(std::max(0.0, mu));
      if (nobs > 1 && std::isfinite(total_var)) {
        prediction_sd.row(g).fill(std::sqrt(std::max(0.0, total_var * (1.0 + 1.0 / nobs))));
        sigma2_out[g] = total_var; method[g] = 1;
        for (int r = 0; r < nobs; ++r) {
          int j = donor_idx[r]; loo(g, j) = (sum1 - y(g, j)) / (nobs - 1.0);
        }
      }
    };
    if (nobs < 2 || !std::isfinite(total_var)) { fallback(); continue; }
    int k = scores.n_cols;
    if (support_adaptive_rank) k = std::min(k, std::max(0, (nobs - 3) / 2));
    if (k < 1 || nobs < std::max(min_target_observed, k + 3) || sse <= 1e-12) {
      fallback(); continue;
    }
    const int p = k + 1;
    arma::mat X(nobs, p, arma::fill::ones); arma::vec target(nobs);
    for (int r = 0; r < nobs; ++r) {
      int j = donor_idx[r]; target[r] = y(g, j);
      for (int a = 0; a < k; ++a) X(r, a + 1) = scores(j, a);
    }
    arma::mat A = X.t() * X;
    for (int a = 1; a < p; ++a) A(a, a) += ridge;
    arma::mat inv;
    if (!dk_inverse(A, inv)) { fallback(); continue; }
    arma::vec beta = inv * X.t() * target;
    arma::vec fitted = X * beta;
    arma::vec hdiag = arma::sum((X * inv) % X, 1);
    arma::vec loo_factor(nobs), loo_null(nobs), d(nobs), centered(nobs);
    for (int r = 0; r < nobs; ++r) {
      double denom = std::max(1e-6, 1.0 - hdiag[r]);
      loo_factor[r] = target[r] - (target[r] - fitted[r]) / denom;
      loo_null[r] = (sum1 - target[r]) / (nobs - 1.0);
      d[r] = loo_factor[r] - loo_null[r];
      centered[r] = target[r] - loo_null[r];
    }
    double den = arma::dot(d, d);
    double q = den > 1e-12 ? arma::dot(d, centered) / den : 0.0;
    q = std::max(0.0, std::min(1.0, q));
    arma::vec loo_cal = loo_null + q * d;
    double loo_bias = arma::mean(loo_cal - target);
    double bias = std::isfinite(bias_kappa) ? nobs / (nobs + bias_kappa) * loo_bias : 0.0;
    if (!std::isfinite(bias)) bias = 0.0;
    loo_cal -= bias;
    double sigma2 = arma::dot(target - loo_cal, target - loo_cal) / nobs;
    if (!std::isfinite(sigma2) || sigma2 < 0) sigma2 = total_var;
    arma::mat Xall(n, p, arma::fill::ones);
    for (int j = 0; j < n; ++j)
      for (int a = 0; a < k; ++a) Xall(j, a + 1) = scores(j, a);
    arma::vec raw = Xall * beta;
    arma::vec hq = arma::sum((Xall * inv) % Xall, 1);
    for (int j = 0; j < n; ++j) {
      prediction(g, j) = std::max(0.0, mu + q * (raw[j] - mu) - bias);
      double lev = (1.0 - q) * (1.0 - q) / nobs + q * q * std::max(0.0, hq[j]);
      prediction_sd(g, j) = std::sqrt(std::max(0.0, sigma2 * (1.0 + lev)));
    }
    for (int r = 0; r < nobs; ++r) loo(g, donor_idx[r]) = loo_cal[r];
    shrinkage[g] = q; bias_calibration[g] = bias; sigma2_out[g] = sigma2;
    factor_rank_used[g] = k; method[g] = q > 0 ? 2 : 1;
  }
  return List::create(
    _["prediction"] = prediction, _["prediction_sd"] = prediction_sd,
    _["loo"] = loo, _["shrinkage"] = shrinkage,
    _["bias_calibration"] = bias_calibration, _["sigma2"] = sigma2_out,
    _["n_observed"] = n_observed, _["factor_rank_used"] = factor_rank_used,
    _["method"] = method
  );
}

static arma::vec dk_simplex_qp(const arma::mat& H, const arma::vec& b,
                               const arma::vec& initial, const double tol) {
  const int k = H.n_rows;
  arma::uvec active = arma::regspace<arma::uvec>(0, k - 1);
  arma::vec w(k, arma::fill::zeros);
  for (int iter = 0; iter < 4 * k + 4; ++iter) {
    arma::mat Ha = H.submat(active, active);
    arma::vec ba = b.elem(active), one(active.n_elem, arma::fill::ones);
    arma::vec Hib, Hi1;
    bool ok1 = arma::solve(Hib, Ha, ba, arma::solve_opts::likely_sympd);
    bool ok2 = arma::solve(Hi1, Ha, one, arma::solve_opts::likely_sympd);
    if (!ok1 || !ok2 || !Hib.is_finite() || !Hi1.is_finite()) break;
    double denom = arma::dot(one, Hi1);
    if (!std::isfinite(denom) || denom <= 0) break;
    double nu = (arma::dot(one, Hib) - 1.0) / denom;
    arma::vec wa = Hib - nu * Hi1;
    arma::uword minpos = wa.index_min(); double minval = wa[minpos];
    if (minval < -tol && active.n_elem > 1) {
      active.shed_row(minpos); continue;
    }
    w.zeros(); w.elem(active) = arma::clamp(wa, 0.0, arma::datum::inf);
    double sw = arma::sum(w); if (sw > 0) w /= sw;
    arma::vec grad = H * w - b;
    double lambda_eq = -arma::mean(grad.elem(active));
    arma::uvec inactive = arma::find(w <= tol);
    if (inactive.n_elem) {
      arma::vec violation = grad.elem(inactive) + lambda_eq;
      arma::uword vp = violation.index_min(); double vv = violation[vp];
      if (vv < -10 * tol) {
        active.insert_rows(active.n_elem, 1); active[active.n_elem - 1] = inactive[vp];
        active = arma::unique(active); continue;
      }
    }
    return w;
  }
  w = arma::clamp(initial, 0.0, arma::datum::inf);
  double sw = arma::sum(w);
  if (sw <= 0) w.fill(1.0 / k); else w /= sw;
  return w;
}

// Exact active-set simplex ridge solver using a precomputed cell Gram matrix.
// [[Rcpp::export]]
arma::mat dk_simplex_weights_cpp(
    const arma::mat& cell_gram,
    const Rcpp::IntegerVector& query_index,
    const Rcpp::IntegerMatrix& candidate_index,
    const arma::mat& geometry_prior,
    const double lambda,
    const double prune_weight,
    const double tol = 1e-10) {
  const int nq = candidate_index.nrow(), K = candidate_index.ncol();
  arma::mat out(nq, K, arma::fill::zeros);
  for (int q = 0; q < nq; ++q) {
    const int c = query_index[q] - 1;
    if (c < 0 || c >= (int)cell_gram.n_rows) continue;
    std::vector<arma::uword> idx; std::vector<int> pos;
    for (int a = 0; a < K; ++a) {
      int j = candidate_index(q, a) - 1;
      if (j >= 0 && j < (int)cell_gram.n_rows && j != c) { idx.push_back(j); pos.push_back(a); }
    }
    if (idx.empty()) continue;
    arma::uvec iu(idx); const int k = idx.size();
    arma::mat H = cell_gram.submat(iu, iu);
    H.diag() += lambda;
    arma::vec prior(k), b(k);
    for (int a = 0; a < k; ++a) prior[a] = std::max(0.0, geometry_prior(q, pos[a]));
    double ps = arma::sum(prior); if (ps <= 0) prior.fill(1.0 / k); else prior /= ps;
    b = cell_gram.submat(iu, arma::uvec({(arma::uword)c})) + lambda * prior;
    arma::vec w = dk_simplex_qp(H, b, prior, tol);
    for (int a = 0; a < k; ++a) if (w[a] < prune_weight) w[a] = 0;
    double sw = arma::sum(w); if (sw <= 0) w = prior; else w /= sw;
    for (int a = 0; a < k; ++a) out(q, pos[a]) = w[a];
  }
  return out;
}

// Cell-wise simplex ridge using only query-observed reliable features.
// A masked target is zero in the query and therefore cannot enter its own
// predictor state.  Each query cell is solved once and all target events in
// that cell share the resulting donor weights.
// [[Rcpp::export]]
Rcpp::List dk_simplex_cell_weights_cpp(
    const arma::mat& features,
    const Rcpp::LogicalMatrix& query_observed,
    const Rcpp::IntegerVector& query_index,
    const Rcpp::IntegerMatrix& candidate_index,
    const arma::mat& geometry_prior,
    const int reliable_features,
    const double lambda,
    const double prune_weight,
    const double tol = 1e-10) {
  const int nq = candidate_index.nrow(), K = candidate_index.ncol();
  if ((int)features.n_rows != query_observed.nrow() ||
      nq != query_observed.ncol() || nq != query_index.size() ||
      nq != (int)geometry_prior.n_rows || K != (int)geometry_prior.n_cols)
    Rcpp::stop("incompatible simplex cell-weight inputs");
  arma::mat out(nq, K, arma::fill::zeros);
  Rcpp::IntegerVector n_features(nq, 0);
  Rcpp::LogicalVector solved(nq, false);
  for (int q = 0; q < nq; ++q) {
    const int c = query_index[q] - 1;
    if (c < 0 || c >= (int)features.n_cols) continue;
    std::vector<arma::uword> rows;
    rows.reserve(std::max(2, reliable_features));
    for (int g = 0; g < query_observed.nrow() &&
                    (int)rows.size() < reliable_features; ++g)
      if (query_observed(g, q) == TRUE) rows.push_back(g);
    if (rows.size() < 2) continue;
    std::vector<arma::uword> idx; std::vector<int> pos;
    for (int a = 0; a < K; ++a) {
      int j = candidate_index(q, a) - 1;
      if (j >= 0 && j < (int)features.n_cols && j != c) {
        idx.push_back(j); pos.push_back(a);
      }
    }
    if (idx.empty()) continue;
    arma::uvec rr(rows), iu(idx); const int k = idx.size();
    arma::mat D = features.submat(rr, iu);
    arma::vec xc = features.submat(rr, arma::uvec({(arma::uword)c}));
    const double scale = 1.0 / rows.size();
    arma::mat H = scale * (D.t() * D);
    H.diag() += lambda;
    arma::vec prior(k), b(k);
    for (int a = 0; a < k; ++a)
      prior[a] = std::max(0.0, geometry_prior(q, pos[a]));
    double ps = arma::sum(prior);
    if (ps <= 0) prior.fill(1.0 / k); else prior /= ps;
    b = scale * (D.t() * xc) + lambda * prior;
    arma::vec w = dk_simplex_qp(H, b, prior, tol);
    for (int a = 0; a < k; ++a) if (w[a] < prune_weight) w[a] = 0;
    double sw = arma::sum(w);
    if (sw <= 0) w = prior; else w /= sw;
    for (int a = 0; a < k; ++a) out(q, pos[a]) = w[a];
    n_features[q] = rows.size(); solved[q] = true;
  }
  return Rcpp::List::create(
    Rcpp::_ ["weights"] = out,
    Rcpp::_ ["n_features"] = n_features,
    Rcpp::_ ["solved"] = solved
  );
}

// Sparse event prediction from one simplex weight vector per query cell.
// [[Rcpp::export]]
Rcpp::List dk_simplex_predict_events_cpp(
    const arma::sp_mat& x,
    const Rcpp::IntegerVector& event_gene,
    const Rcpp::IntegerVector& event_cell,
    const Rcpp::IntegerMatrix& candidate_index,
    const arma::mat& cell_weights,
    const double positive_mass_tol = 1e-8) {
  const int E = event_gene.size(), K = candidate_index.ncol();
  if (event_cell.size() != E || candidate_index.nrow() != (int)x.n_cols ||
      cell_weights.n_rows != x.n_cols || cell_weights.n_cols != (arma::uword)K)
    Rcpp::stop("incompatible simplex event-prediction inputs");
  Rcpp::NumericVector prediction(E, NA_REAL), prediction_sd(E, NA_REAL),
    sum_weight_error(E, NA_REAL), min_weight(E, NA_REAL);
  for (int e = 0; e < E; ++e) {
    const int g = event_gene[e] - 1, c = event_cell[e] - 1;
    if (g < 0 || g >= (int)x.n_rows || c < 0 || c >= (int)x.n_cols) continue;
    double mass = 0.0, num = 0.0, sw = 0.0,
      minw = std::numeric_limits<double>::infinity();
    for (int a = 0; a < K; ++a) {
      const int d = candidate_index(c, a) - 1; const double w = cell_weights(c, a);
      if (d < 0 || d >= (int)x.n_cols || !std::isfinite(w) || w <= 0) continue;
      sw += w; minw = std::min(minw, w);
      const double val = x(g, d);
      if (std::isfinite(val) && val > 0) { mass += w; num += w * val; }
    }
    if (!(mass > positive_mass_tol) || !(sw > 0)) continue;
    const double mu = num / mass;
    double vnum = 0.0, wp2 = 0.0;
    for (int a = 0; a < K; ++a) {
      const int d = candidate_index(c, a) - 1; const double w = cell_weights(c, a);
      if (d < 0 || d >= (int)x.n_cols || !std::isfinite(w) || w <= 0) continue;
      const double val = x(g, d);
      if (std::isfinite(val) && val > 0) {
        const double wp = w / mass;
        vnum += wp * (val - mu) * (val - mu); wp2 += wp * wp;
      }
    }
    const double denom = 1.0 - wp2;
    const double vv = denom > positive_mass_tol ? vnum / denom : 0.0;
    prediction[e] = mu;
    prediction_sd[e] = std::sqrt(std::max(0.0, vv * (1.0 + wp2)));
    sum_weight_error[e] = std::abs(sw - 1.0);
    min_weight[e] = std::isfinite(minw) ? minw : NA_REAL;
  }
  return Rcpp::List::create(
    Rcpp::_ ["prediction"] = prediction,
    Rcpp::_ ["prediction_sd"] = prediction_sd,
    Rcpp::_ ["sum_weight_error"] = sum_weight_error,
    Rcpp::_ ["min_weight"] = min_weight
  );
}
