// Brute-force max-norm nearest-neighbour kernels for Rcausality.
// All distances are Chebyshev (max-norm) distances between rows; matrices
// are copied once into row-major buffers so that the O(n^2 d) loops are
// cache-friendly.  Ties and the "strictly closer than the k-th neighbour"
// convention follow the pure-R reference implementations in R/knn_cmi.R.
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <functional>
using namespace Rcpp;

namespace {

struct RowMat {
  int n, d;
  std::vector<double> v;           // row-major n x d
  RowMat(const NumericMatrix& M) : n(M.nrow()), d(M.ncol()), v((size_t)M.nrow() * M.ncol()) {
    for (int i = 0; i < n; ++i) for (int j = 0; j < d; ++j) v[(size_t)i * d + j] = M(i, j);
  }
  inline double dist(int i, int j) const {
    double m = 0.0;
    const double* a = &v[(size_t)i * d]; const double* b = &v[(size_t)j * d];
    for (int k = 0; k < d; ++k) { double t = std::fabs(a[k] - b[k]); if (t > m) m = t; }
    return m;
  }
};

// k-th smallest of the distances from i to all j != i (k >= 1)
inline double kth_row(const std::vector<double>& row, int i, int k, std::vector<double>& buf) {
  buf.clear();
  for (int j = 0; j < (int)row.size(); ++j) if (j != i) buf.push_back(row[j]);
  if (k > (int)buf.size()) k = (int)buf.size();
  std::nth_element(buf.begin(), buf.begin() + (k - 1), buf.end());
  return buf[k - 1];
}

} // namespace

// Frenzel-Pompe / KSG estimator of I(X; Y | Z) in nats (Z may have 0 columns).
// [[Rcpp::export]]
double cmi_knn_cpp(NumericMatrix X, NumericMatrix Y, NumericMatrix Z, int k) {
  RowMat x(X), y(Y), z(Z);
  const int n = x.n; const bool hasz = z.d > 0;
  std::vector<double> dx(n), dy(n), dz(n), dj(n), buf;
  double acc = 0.0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      dx[j] = x.dist(i, j); dy[j] = y.dist(i, j); dz[j] = hasz ? z.dist(i, j) : 0.0;
      dj[j] = std::max(dx[j], std::max(dy[j], dz[j]));
    }
    double eps = kth_row(dj, i, k, buf);
    int nxz = -1, nyz = -1, nz = -1;             // the point itself is subtracted
    for (int j = 0; j < n; ++j) {
      if (std::max(dx[j], dz[j]) < eps) ++nxz;
      if (std::max(dy[j], dz[j]) < eps) ++nyz;
      if (dz[j] < eps) ++nz;
    }
    acc += R::digamma(nxz + 1.0) + R::digamma(nyz + 1.0) - R::digamma(nz + 1.0);
  }
  return R::digamma((double)k) - acc / n;
}

// indices (1-based) of the k nearest neighbours of every row (self excluded),
// ordered by increasing max-norm distance
// [[Rcpp::export]]
IntegerMatrix knn_idx_cpp(NumericMatrix M, int k) {
  RowMat m(M); const int n = m.n;
  if (k > n - 1) k = n - 1;
  IntegerMatrix out(n, k);
  std::vector<std::pair<double, int> > row(n);
  for (int i = 0; i < n; ++i) {
    int cnt = 0;
    for (int j = 0; j < n; ++j) if (j != i) row[cnt++] = std::make_pair(m.dist(i, j), j);
    std::partial_sort(row.begin(), row.begin() + k, row.begin() + cnt);
    for (int l = 0; l < k; ++l) out(i, l) = row[l].second + 1;
  }
  return out;
}

// Conditional-variance kernels of the df method.  C: constraint matrix
// (n x d, d may be 0), P: present values (n x m).  For every centre i the
// samples whose constraint vector lies within the window are pooled.
// Returns the mean over the columns of P of log(pooled within-window variance)
// with adaptive windows containing the k nearest neighbours (centre included,
// unbiased within-group variance).
// [[Rcpp::export]]
double cond_var_knn_cpp(NumericMatrix C, NumericMatrix P, int k) {
  RowMat c(C), p(P); const int n = c.n, m = p.d;
  std::vector<double> num(m, 0.0); double den = 0.0;
  if (c.d == 0) {                                  // no constraint: plain variance
    double out = 0.0;
    for (int j = 0; j < m; ++j) {
      double s = 0, s2 = 0;
      for (int i = 0; i < n; ++i) { double v = p.v[(size_t)i * m + j]; s += v; s2 += v * v; }
      double var = (s2 - s * s / n) / (n - 1.0);
      out += std::log(std::max(var, 1e-300));
    }
    return out / m;
  }
  std::vector<double> row(n), buf;
  std::vector<int> idx; idx.reserve(n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) row[j] = c.dist(i, j);
    double eps = kth_row(row, i, k, buf);
    idx.clear();
    for (int j = 0; j < n; ++j) if (row[j] <= eps) idx.push_back(j);
    const int cnt = (int)idx.size();
    if (cnt < 2) continue;
    for (int j = 0; j < m; ++j) {
      double s = 0, s2 = 0;
      for (int l = 0; l < cnt; ++l) { double v = p.v[(size_t)idx[l] * m + j]; s += v; s2 += v * v; }
      double ss = s2 - s * s / cnt;               // (cnt - 1) * unbiased variance
      num[j] += ss;
    }
    den += cnt - 1.0;
  }
  double out = 0.0;
  for (int j = 0; j < m; ++j) out += std::log(std::max(num[j] / den, 1e-300));
  return out / m;
}

// Fixed-width windows (radius r): pooled variance over the centres with at
// least min_n samples; returns c(sigma2, mean window count) or NA when fewer
// than 5 centres qualify.
// [[Rcpp::export]]
NumericVector cond_var_radius_cpp(NumericMatrix C, NumericMatrix P, double r, int min_n) {
  RowMat c(C), p(P); const int n = c.n, m = p.d;
  std::vector<double> num(m, 0.0); double den = 0.0, cntsum = 0.0; int ok = 0;
  std::vector<int> idx; idx.reserve(n);
  for (int i = 0; i < n; ++i) {
    idx.clear();
    for (int j = 0; j < n; ++j) if (c.dist(i, j) < r) idx.push_back(j);
    const int cnt = (int)idx.size();
    if (cnt < min_n) continue;
    ++ok; cntsum += cnt;
    for (int j = 0; j < m; ++j) {
      double s = 0, s2 = 0;
      for (int l = 0; l < cnt; ++l) { double v = p.v[(size_t)idx[l] * m + j]; s += v; s2 += v * v; }
      num[j] += s2 - s * s / cnt;
    }
    den += cnt - 1.0;
  }
  NumericVector out(2);
  if (ok < 5) { out[0] = NA_REAL; out[1] = 0; return out; }
  double v = 0; for (int j = 0; j < m; ++j) v += num[j] / den;
  out[0] = v / m; out[1] = cntsum / ok;
  return out;
}

// Local permutation (Runge 2018): every index is mapped to one of its kp
// nearest neighbours (1-based matrix nn, rows ordered by distance), each
// target used at most once when possible.  `order` is the (1-based) random
// visiting order, `u` uniform draws for the choices, both supplied from R so
// that the RNG state stays under R's control.
// [[Rcpp::export]]
IntegerVector local_perm_cpp(IntegerMatrix nn, int kp, IntegerVector order, NumericVector u) {
  const int n = nn.nrow();
  IntegerVector perm(n); std::vector<char> used(n, 0); std::vector<int> cand; cand.reserve(kp);
  for (int t = 0; t < n; ++t) {
    int i = order[t] - 1;
    cand.clear();
    for (int l = 0; l < kp; ++l) { int j = nn(i, l) - 1; if (!used[j]) cand.push_back(j); }
    int j;
    if (!cand.empty()) j = cand[(int)(u[t] * cand.size()) % cand.size()];
    else j = nn(i, (int)(u[t] * kp) % kp) - 1;
    perm[i] = j + 1; used[j] = 1;
  }
  return perm;
}

// p_markov kernel of the Markov test (Algorithm 2 of Stippinger et al. 2023).
// S: state matrix (N x d, standardised); idx: ni x K matrix of 1-based time
// indices of the neighbourhood of every interior centre; fut: 1-based row of
// the "future" of every time index (t + 1 for the original series, the future
// of a locally permuted neighbour for the surrogates); Up, Uf: ni x d random
// projection directions for pasts and futures; b: number of bins.
// For every centre the pasts S[idx - 1] and futures S[fut[idx]] are projected,
// binned into b quantile bins (ties broken by the supplied uniform jitter
// `jit`, length ni*K), cross-tabulated, Yates-corrected and tested with the
// modified log-likelihood statistic.  Returns the ni local p-values.
// [[Rcpp::export]]
NumericVector markov_local_p_cpp(NumericMatrix S, IntegerMatrix idx, IntegerVector fut, NumericMatrix Up,
                                 NumericMatrix Uf, int b, NumericVector jit) {
  const int ni = idx.nrow(), K = idx.ncol(), d = S.ncol();
  NumericVector out(ni);
  std::vector<std::pair<double,int> > pp(K), ff(K);
  std::vector<int> bp(K), bf(K);
  std::vector<double> tab((size_t)b * b), rs(b), cs(b);
  for (int i = 0; i < ni; ++i) {
    for (int l = 0; l < K; ++l) {
      int t = idx(i, l) - 1; int tf = fut[t] - 1;
      double sp = 0, sf = 0;
      for (int j = 0; j < d; ++j) { sp += S(t - 1, j) * Up(i, j); sf += S(tf, j) * Uf(i, j); }
      double e = jit[(size_t)i * K + l] * 1e-9;
      pp[l] = std::make_pair(sp + e, l); ff[l] = std::make_pair(sf + e, l);
    }
    std::sort(pp.begin(), pp.end()); std::sort(ff.begin(), ff.end());
    for (int r = 0; r < K; ++r) {              // quantile bins: ceiling((r+1) * b / K)
      int bin = (int)std::ceil((r + 1.0) * b / K) - 1; if (bin >= b) bin = b - 1;
      bp[pp[r].second] = bin; bf[ff[r].second] = bin;
    }
    std::fill(tab.begin(), tab.end(), 0.0); std::fill(rs.begin(), rs.end(), 0.0); std::fill(cs.begin(), cs.end(), 0.0);
    for (int l = 0; l < K; ++l) { tab[bp[l] * b + bf[l]] += 1.0; rs[bp[l]] += 1.0; cs[bf[l]] += 1.0; }
    int nr = 0, nc = 0; for (int a = 0; a < b; ++a) { if (rs[a] > 0) ++nr; if (cs[a] > 0) ++nc; }
    double stat = 0.0;
    for (int a = 0; a < b; ++a) for (int c = 0; c < b; ++c) {
      double E = rs[a] * cs[c] / K; if (E <= 0) continue;
      double O = tab[a * b + c], dlt = E - O;
      double Oy = O + (dlt > 0 ? 1.0 : (dlt < 0 ? -1.0 : 0.0)) * std::min(0.5, std::fabs(dlt));
      stat += 2.0 * E * std::log(E / std::max(Oy, 1e-12));
    }
    int df = (nr - 1) * (nc - 1);
    double p = df > 0 ? R::pchisq(std::max(stat, 0.0), (double)df, 0, 0) : 1.0;
    out[i] = std::max(p, 1e-300);
  }
  return out;
}

// sorted Euclidean distances to the kmax nearest neighbours of every row
// (self excluded); used by the intrinsic-dimension estimator
// [[Rcpp::export]]
NumericMatrix knn_dist_cpp(NumericMatrix M, int kmax) {
  RowMat m(M); const int n = m.n, d = m.d;
  if (kmax > n - 1) kmax = n - 1;
  NumericMatrix out(n, kmax);
  std::vector<double> row(n);
  for (int i = 0; i < n; ++i) {
    int cnt = 0;
    const double* a = &m.v[(size_t)i * d];
    for (int j = 0; j < n; ++j) {
      if (j == i) continue;
      const double* b = &m.v[(size_t)j * d]; double s = 0.0;
      for (int l = 0; l < d; ++l) { double t = a[l] - b[l]; s += t * t; }
      row[cnt++] = std::sqrt(s);
    }
    std::partial_sort(row.begin(), row.begin() + kmax, row.begin() + cnt);
    for (int l = 0; l < kmax; ++l) out(i, l) = row[l];
  }
  return out;
}

// Non-sequential recursive pair substitution (NSRPS) compression rates
// (CCDH, Stippinger et al.): at every step the most frequent non-overlapping
// consecutive pair (ties: the one occurring first) is replaced by a new
// symbol; returns |G^k x| / |x| for k = 1..K.
// [[Rcpp::export]]
NumericVector nsrps_cr_cpp(IntegerVector x, int K) {
  std::vector<long long> s(x.begin(), x.end());
  const double n0 = (double)s.size();
  long long next_sym = 0; for (size_t i = 0; i < s.size(); ++i) next_sym = std::max(next_sym, s[i]);
  NumericVector out(K);
  for (int k = 0; k < K; ++k) {
    const size_t n = s.size();
    if (n < 2) { out[k] = n / n0; continue; }
    // check whether all symbols are equal (entropy zero): stop
    bool allsame = true; for (size_t i = 1; i < n; ++i) if (s[i] != s[0]) { allsame = false; break; }
    if (allsame) { out[k] = n / n0; continue; }
    // count non-overlapping occurrences of every pair
    std::vector<std::pair<long long,long long> > keys; keys.reserve(n);
    std::vector<int> cnt; std::vector<size_t> first;
    // map pair -> index via sorting-free approach: use std::vector + linear map (symbol space small)
    // use an unordered map keyed on combined value
    struct H { size_t operator()(const std::pair<long long,long long>& p) const {
      return std::hash<long long>()(p.first * 1000003LL + p.second); } };
    std::vector<size_t> last_end;
    std::unordered_map<std::pair<long long,long long>, size_t, H> idx;
    for (size_t i = 0; i + 1 < n; ++i) {
      std::pair<long long,long long> p(s[i], s[i + 1]);
      auto it = idx.find(p);
      if (it == idx.end()) {
        idx.emplace(p, keys.size()); keys.push_back(p); cnt.push_back(1); first.push_back(i); last_end.push_back(i + 1);
      } else {
        size_t j = it->second;
        if (last_end[j] >= i) continue;           // overlaps the previous counted occurrence
        cnt[j]++; last_end[j] = i + 1;
      }
    }
    size_t best = 0;
    for (size_t j = 1; j < keys.size(); ++j)
      if (cnt[j] > cnt[best] || (cnt[j] == cnt[best] && first[j] < first[best])) best = j;
    const long long a = keys[best].first, b = keys[best].second;
    ++next_sym;
    std::vector<long long> t; t.reserve(n);
    for (size_t i = 0; i < n; ) {
      if (i + 1 < n && s[i] == a && s[i + 1] == b) { t.push_back(next_sym); i += 2; }
      else { t.push_back(s[i]); ++i; }
    }
    s.swap(t);
    out[k] = s.size() / n0;
  }
  return out;
}

// ---- df-entropy method (Jakovac et al.): binary conditional entropies -------
// Conditional entropy (bits) of `target` (integer states 0..S-1) given the
// pattern formed by the 0/1 columns of `hist` (rows aligned with target);
// zero-frequency states are omitted (Appendix F.0.2).
// [[Rcpp::export]]
double cond_entropy_bits_cpp(IntegerVector target, IntegerMatrix hist, int n_states) {
  const int n = target.size(), d = hist.ncol();
  if (hist.nrow() != n) Rcpp::stop("target and history must have the same number of rows");
  if (d > 24) Rcpp::stop("too many history bits (at most 24)");
  if (n_states < 1) Rcpp::stop("n_states must be positive");
  const long long npat = 1LL << d;
  std::vector<double> cnt((size_t)npat * n_states, 0.0), tot((size_t)npat, 0.0);
  for (int i = 0; i < n; ++i) {
    long long code = 0;
    if (target[i] == NA_INTEGER || target[i] < 0 || target[i] >= n_states) Rcpp::stop("target values must be in 0..n_states-1 without NA");
    for (int j = 0; j < d; ++j) { if (hist(i, j) == NA_INTEGER) Rcpp::stop("NA in the history"); if (hist(i, j)) code |= (1LL << j); }
    cnt[(size_t)code * n_states + target[i]] += 1.0; tot[(size_t)code] += 1.0;
  }
  double H = 0.0;
  for (long long c = 0; c < npat; ++c) {
    if (tot[c] <= 0) continue;
    for (int s = 0; s < n_states; ++s) {
      double p = cnt[(size_t)c * n_states + s];
      if (p > 0) H -= p / n * std::log2(p / tot[c]);
    }
  }
  return H;
}

// Coupled binary Markov chains (training data of the df-entropy method).
// x[t+1] ~ Bernoulli(px[code]) where code packs (x[t], ..., x[t-mx+1]) and,
// when coupled, the current state of the driver(s) (y[t] and/or z[t]);
// analogously for y; z is an autonomous chain of memory mz (hidden common
// cause). `u` supplies the uniform draws (3 * n) from R's RNG.
// [[Rcpp::export]]
IntegerMatrix sim_coupled_binary_cpp(int n, int mx, int my, int mz, bool y_to_x, bool x_to_y,
                                     bool z_to_x, bool z_to_y, NumericVector px, NumericVector py,
                                     NumericVector pz, NumericVector u) {
  IntegerMatrix out(n, 3);
  std::vector<int> x(n), y(n), z(n);
  int m = std::max(mx, std::max(my, mz));
  for (int t = 0; t < m; ++t) { x[t] = u[3 * t] < 0.5; y[t] = u[3 * t + 1] < 0.5; z[t] = u[3 * t + 2] < 0.5; }
  auto code = [](const std::vector<int>& v, int t, int mm) {
    long long c = 0;
    for (int j = 0; j < mm; ++j) { if (v[t - j]) c |= (1LL << j); }
    return c;
  };
  for (int t = m - 1; t + 1 < n; ++t) {
    long long cz = code(z, t, mz);
    z[t + 1] = u[3 * (t + 1) + 2] < pz[cz];
    long long cx = code(x, t, mx); int bit = mx;
    if (y_to_x) { if (y[t]) cx |= (1LL << bit); ++bit; }
    if (z_to_x) { if (z[t]) cx |= (1LL << bit); ++bit; }
    x[t + 1] = u[3 * (t + 1)] < px[cx];
    long long cy = code(y, t, my); bit = my;
    if (x_to_y) { if (x[t]) cy |= (1LL << bit); ++bit; }
    if (z_to_y) { if (z[t]) cy |= (1LL << bit); ++bit; }
    y[t + 1] = u[3 * (t + 1) + 1] < py[cy];
  }
  for (int t = 0; t < n; ++t) { out(t, 0) = x[t]; out(t, 1) = y[t]; out(t, 2) = z[t]; }
  return out;
}

// ---- transfer entropy / CCM helpers ----------------------------------------
// Conditional entropy (bits) of `target` (states 0..n_states-1) given the
// pattern of the integer columns of `hist`, each in 0..base-1.
// [[Rcpp::export]]
double cond_entropy_sym_cpp(IntegerVector target, IntegerMatrix hist, int n_states, int base) {
  const int n = target.size(), d = hist.ncol();
  if (hist.nrow() != n) Rcpp::stop("target and history must have the same number of rows");
  std::unordered_map<long long, std::vector<double> > tab;
  for (int i = 0; i < n; ++i) {
    long long code = 0, mult = 1;
    if (target[i] == NA_INTEGER || target[i] < 0 || target[i] >= n_states) Rcpp::stop("target values must be in 0..n_states-1 without NA");
    for (int j = 0; j < d; ++j) { if (hist(i, j) == NA_INTEGER || hist(i, j) < 0) Rcpp::stop("NA or negative symbol in the history"); code += (long long)hist(i, j) * mult; mult *= base; }
    std::vector<double>& v = tab[code];
    if (v.empty()) v.assign(n_states + 1, 0.0);
    v[target[i]] += 1.0; v[n_states] += 1.0;
  }
  double H = 0.0;
  for (auto& kv : tab) {
    const std::vector<double>& v = kv.second; double tot = v[n_states];
    for (int s = 0; s < n_states; ++s) if (v[s] > 0) H -= v[s] / n * std::log2(v[s] / tot);
  }
  return H;
}

// Simplex projection (Sugihara & May 1990) used by CCM: for every row of
// `query` find its k nearest library points (Euclidean, rows `lib` of M,
// excluding the query's own row when `self` gives it), and return the
// exponentially weighted average of `value` over those neighbours.
// [[Rcpp::export]]
NumericVector simplex_predict_cpp(NumericMatrix M, IntegerVector lib, IntegerVector query,
                                  NumericVector value, int k) {
  RowMat m(M); const int nl = lib.size(), nq = query.size();
  NumericVector out(nq);
  std::vector<std::pair<double,int> > row; row.reserve(nl);
  for (int q = 0; q < nq; ++q) {
    const int i = query[q] - 1;
    row.clear();
    for (int l = 0; l < nl; ++l) { int j = lib[l] - 1; if (j == i) continue;
      double s = 0; for (int c = 0; c < m.d; ++c) { double t = m.v[(size_t)i * m.d + c] - m.v[(size_t)j * m.d + c]; s += t * t; }
      row.push_back(std::make_pair(std::sqrt(s), j)); }
    int kk = std::min(k, (int)row.size());
    if (kk == 0) { out[q] = NA_REAL; continue; }
    std::partial_sort(row.begin(), row.begin() + kk, row.end());
    double d1 = std::max(row[0].first, 1e-12), wsum = 0, acc = 0;
    for (int l = 0; l < kk; ++l) { double w = std::exp(-row[l].first / d1); wsum += w; acc += w * value[row[l].second]; }
    out[q] = acc / wsum;
  }
  return out;
}

// ---- recurrence plots (Hirata & Aihara 2010) and correlation sums -----------
// Pairwise Euclidean distances of the rows of M (upper triangle, row-major),
// used for recurrence thresholds and correlation sums.
// [[Rcpp::export]]
NumericVector pair_dist_cpp(NumericMatrix M) {
  RowMat m(M); const int n = m.n; const long long np = (long long)n * (n - 1) / 2;
  NumericVector out(np); long long c = 0;
  for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) {
    double s = 0; for (int l = 0; l < m.d; ++l) { double t = m.v[(size_t)i * m.d + l] - m.v[(size_t)j * m.d + l]; s += t * t; }
    out[c++] = std::sqrt(s);
  }
  return out;
}

// Recurrence-plot statistics for two embedded series with a loose and a
// tight threshold per series (recurrence = distance <= eps), restricted to
// pairs with |i - j| >= w (Theiler window); for surrogates the second series
// is read at the circularly shifted pair (i + s, j + s).  Returns
// c(n_pairs, rx_loose, ry_loose, rx_tight, ry_tight, joint_loose,
//   n[R_X^loose & R_Y^tight], n[R_Y^loose & R_X^tight]).
// [[Rcpp::export]]
NumericVector recurrence_counts_cpp(NumericMatrix X, NumericMatrix Y, double ex_loose, double ex_tight,
                                    double ey_loose, double ey_tight, int w, int s) {
  RowMat x(X), y(Y); const int n = x.n;
  double np = 0, rxl = 0, ryl = 0, rxt = 0, ryt = 0, rj = 0, cxy = 0, cyx = 0;
  if (w < 1) w = 1;                                   // the diagonal i == j is never a recurrence
  for (int i = 0; i < n; ++i) for (int j = i + w; j < n; ++j) {
    double dx = 0; for (int l = 0; l < x.d; ++l) { double t = x.v[(size_t)i * x.d + l] - x.v[(size_t)j * x.d + l]; dx += t * t; }
    int is = (i + s) % n, js = (j + s) % n;
    double dy = 0; for (int l = 0; l < y.d; ++l) { double t = y.v[(size_t)is * y.d + l] - y.v[(size_t)js * y.d + l]; dy += t * t; }
    dx = std::sqrt(dx); dy = std::sqrt(dy);
    bool axl = dx <= ex_loose, axt = dx <= ex_tight, ayl = dy <= ey_loose, ayt = dy <= ey_tight;
    np += 1; rxl += axl; ryl += ayl; rxt += axt; ryt += ayt; rj += (axl && ayl);
    cxy += (axl && ayt); cyx += (ayl && axt);
  }
  return NumericVector::create(np, rxl, ryl, rxt, ryt, rj, cxy, cyx);
}
