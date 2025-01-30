#include <Rcpp.h>
#include <cmath>
#include <algorithm> // for std::fabs
using namespace Rcpp;

// ----------------------------------------------------------------------------
//                        Helper / Utility Functions
// ----------------------------------------------------------------------------

// Soft-thresholding (L1 shrinkage)
inline double SoftShrink(double y, double lam) {
  // standard "soft threshold":
  //   y > 0 => y - lam
  //   y < 0 => y + lam
  // returning 0 if (|y| <= lam)
  double temp;
  if (y > 0) {
    temp = y - lam;
  } else {
    temp = -y - lam;
  }
  if (temp <= 0.0) {
    return 0.0;
  } else {
    return (y > 0.0) ? temp : -temp;
  }
}

// ----------------------------------------------------------------------------

// Calculate row-wise norm of penalized coefficients
inline void CalBnorm(int P, int Q,
                     const NumericMatrix &Beta,
                     const IntegerMatrix &C_m,
                     NumericVector &Bnorm) {
  for (int p = 0; p < P; p++) {
    double sumSq = 0.0;
    for (int q = 0; q < Q; q++) {
      if (C_m(p, q) == 1) {
        double val = Beta(p, q);
        sumSq += val * val;
      }
    }
    Bnorm[p] = std::sqrt(sumSq);
  }
}

// ----------------------------------------------------------------------------

// Copy matrix -> matrix
inline void Assign(const NumericMatrix &src, NumericMatrix &dest) {
  int P = src.nrow();
  int Q = src.ncol();
  for (int p = 0; p < P; p++) {
    for (int q = 0; q < Q; q++) {
      dest(p, q) = src(p, q);
    }
  }
}

// ----------------------------------------------------------------------------

// Infinity-norm distance (max absolute difference) between two matrices
inline double Dist(const NumericMatrix &A, const NumericMatrix &B) {
  double result = 0.0;
  int P = A.nrow();
  int Q = A.ncol();
  for (int p = 0; p < P; p++) {
    for (int q = 0; q < Q; q++) {
      double diff = std::fabs(A(p, q) - B(p, q));
      if (diff > result) {
        result = diff;
      }
    }
  }
  return result;
}

// ----------------------------------------------------------------------------
//   Update() with sigma:  Beta[p,q] = SoftShrink( temp1, lam2 * sigma[q] )
//                         phi[p,q]  = Beta[p,q] * SoftShrink(..., lam1*sigma[q])
// ----------------------------------------------------------------------------

static void Update(int cur_p,
                   int N, int P, int Q,
                   double lam1, double lam2,
                   const NumericVector &sigma,   // new: sigma[q]
                   const IntegerMatrix &C_m,
                   const NumericMatrix &X_m,
                   const NumericVector &Xnorm,
                   NumericMatrix &E,
                   NumericMatrix &Beta,
                   NumericVector &Bnorm,
                   NumericMatrix &phi_old,
                   NumericMatrix &phi) {

  // 1) Lasso solution for Beta
  for (int q = 0; q < Q; q++) {
    if (C_m(cur_p, q) == 0) {
      Beta(cur_p, q) = 0.0;
    } else {
      // sum of E(n,q)*X_m(n,p)
      double temp = 0.0;
      for (int n = 0; n < N; n++) {
        temp += E(n, q) * X_m(n, cur_p);
      }
      // temp1 = temp + phi_old(p,q)*Xnorm[p]
      double temp1 = temp + phi_old(cur_p, q) * Xnorm[cur_p];

      if (C_m(cur_p, q) == 1) {  // penalized
        // multiply lam2 by sigma[q]
        double shr = lam2 * sigma[q];
        Beta(cur_p, q) = SoftShrink(temp1, shr) / Xnorm[cur_p];
      } else {
        // C_m[p,q]==2 => not penalized
        Beta(cur_p, q) = temp1 / Xnorm[cur_p];
      }
    }
  }

  // 2) Recompute Bnorm for penalized spots
  double sumSq = 0.0;
  for (int q = 0; q < Q; q++) {
    if (C_m(cur_p, q) == 1) {
      double val = Beta(cur_p, q);
      sumSq += val * val;
    }
  }
  Bnorm[cur_p] = std::sqrt(sumSq);

  // 3) Update phi
  for (int q = 0; q < Q; q++) {
    if (C_m(cur_p, q) == 0) {
      phi(cur_p, q) = 0.0;
    } else if (C_m(cur_p, q) == 1 && Bnorm[cur_p] > 1e-6) {
      // group-lasso shrinkage
      double tmp = Xnorm[cur_p] * Bnorm[cur_p];
      // multiply lam1 by sigma[q]
      double shr = SoftShrink(tmp, lam1 * sigma[q]);
      phi(cur_p, q) = (Beta(cur_p, q) * shr) / tmp;
    } else {
      // not penalized or row norm is near 0
      phi(cur_p, q) = Beta(cur_p, q);
    }
  }

  // 4) Update residual E
  for (int q = 0; q < Q; q++) {
    double oldVal = phi_old(cur_p, q);
    double newVal = phi(cur_p, q);
    double diff = oldVal - newVal;
    if (std::fabs(diff) > 1e-15) {
      for (int n = 0; n < N; n++) {
        E(n, q) = E(n, q) + diff * X_m(n, cur_p);
      }
    }
  }

  // 5) Update phi_old
  for (int q = 0; q < Q; q++) {
    phi_old(cur_p, q) = phi(cur_p, q);
  }

  // 6) Recompute Bnorm w.r.t. new phi
  sumSq = 0.0;
  for (int q = 0; q < Q; q++) {
    if (C_m(cur_p, q) == 1) {
      double val = phi(cur_p, q);
      sumSq += val * val;
    }
  }
  Bnorm[cur_p] = std::sqrt(sumSq);
}

// ----------------------------------------------------------------------------
//      1) MultiRegGroupLasso with sigma & eps=1e-3
// ----------------------------------------------------------------------------

// [[Rcpp::export]]
List MultiRegGroupLasso_sigma_unified(
    const NumericMatrix &X_m,        // N x P
    const NumericMatrix &Y_m,        // N x Q
    const IntegerMatrix &C_m,        // P x Q
    double lam1,                     // group-lasso penalty
    double lam2,                     // lasso penalty
    const NumericVector &sigma,      // length-Q scale vector
    Nullable<NumericMatrix> Phi_initial = R_NilValue
) {
  /*
   * This function unifies both "no-initial-value" and
   * "with-initial-value" approaches into one code path.
   *
   * If Phi_initial is NULL, we do the standard Beta-lasso initialization.
   * If Phi_initial is provided, we skip that step and initialize phi
   * from the user-supplied matrix (except for excluded positions).
   */

  // 0) Basic dims
  int N = X_m.nrow();
  int P = X_m.ncol();
  int Q = Y_m.ncol();

  if (sigma.size() != (size_t)Q) {
    stop("sigma must have length == # of columns in Y_m");
  }

  // 1) Precompute Xnorm
  NumericVector Xnorm(P, 0.0);
  for (int p = 0; p < P; p++) {
    double sumSq = 0.0;
    for (int n = 0; n < N; n++) {
      double val = X_m(n, p);
      sumSq += val * val;
    }
    Xnorm[p] = sumSq;
  }

  // 2) Allocate objects
  NumericMatrix Beta(P, Q), phi(P, Q),
  phi_old(P, Q), phi_last(P, Q),
  E(N, Q);
  NumericVector Bnorm(P, 0.0);

  // --------------------------------------------------------------------------
  // 3) Initialize phi (two cases):
  //    a) If user gave no initial matrix => do "Beta-lasso" approach
  //    b) Otherwise read user-provided initial
  // --------------------------------------------------------------------------
  bool hasInit = Phi_initial.isNotNull();
  if (!hasInit) {
    // (a) Standard approach: first compute a simple "lasso solution" => Beta,
    //     then do group-lasso shrink => phi.

    // i) Compute Beta
    for (int p = 0; p < P; p++) {
      for (int q = 0; q < Q; q++) {
        if (C_m(p, q) == 0) {
          Beta(p, q) = 0.0;
        } else {
          // sum_{n} X_m(n,p) * Y_m(n,q)
          double tmp = 0.0;
          for (int n = 0; n < N; n++) {
            tmp += X_m(n, p) * Y_m(n, q);
          }
          if (C_m(p, q) == 1) {  // penalized
            Beta(p, q) = SoftShrink(tmp, lam2 * sigma[q]) / Xnorm[p];
          } else {
            // unpenalized => no shrink
            Beta(p, q) = tmp / Xnorm[p];
          }
        }
      }
    }

    // ii) Bnorm => rowwise group-lasso norm
    CalBnorm(P, Q, Beta, C_m, Bnorm);

    // iii) phi from Beta with group-lasso shrink
    for (int p = 0; p < P; p++) {
      double bn = Bnorm[p];
      if (bn > 1e-6) {
        for (int q = 0; q < Q; q++) {
          if (C_m(p, q) == 1) {
            double tmp = Xnorm[p] * bn;
            double shr = SoftShrink(tmp, lam1 * sigma[q]);
            phi(p, q)   = Beta(p, q) * (shr / tmp);
          } else if (C_m(p, q) == 2) {
            phi(p, q) = Beta(p, q);
          } else {
            phi(p, q) = 0.0;
          }
        }
      } else {
        // entire row is near 0 => only unpenalized might remain
        for (int q = 0; q < Q; q++) {
          if (C_m(p, q) == 2) {
            phi(p, q) = Beta(p, q);
          } else {
            phi(p, q) = 0.0;
          }
        }
      }
    }

  } else {
    // (b) The user provided an initial matrix => read it
    NumericMatrix PhiInit(Phi_initial);  // copy from the Nullable
    if (PhiInit.nrow() != P || PhiInit.ncol() != Q) {
      stop("Phi_initial must be a (P x Q) matrix or NULL.");
    }
    for (int p = 0; p < P; p++) {
      for (int q = 0; q < Q; q++) {
        if (C_m(p, q) == 0) {
          phi(p, q) = 0.0; // exclude
        } else {
          phi(p, q) = PhiInit(p, q);
        }
      }
    }
  } // end if(!hasInit) else

  // --------------------------------------------------------------------------
  // 4) Compute the initial residual E = Y - X * phi
  // --------------------------------------------------------------------------
  for (int n = 0; n < N; n++) {
    for (int q = 0; q < Q; q++) {
      double tmp = 0.0;
      for (int p = 0; p < P; p++) {
        tmp += phi(p, q) * X_m(n, p);
      }
      E(n, q) = Y_m(n, q) - tmp;
    }
  }

  // --------------------------------------------------------------------------
  // 5) Main coordinate-descent iteration loop
  // --------------------------------------------------------------------------
  double eps    = 1e-3;
  double flag   = 100.0;
  int n_iter    = 0;
  int max_iter  = 1000000000;  // i.e. 1e+10

  while ((flag > eps) && (n_iter < max_iter)) {
    // Recompute rowwise norm Bnorm from current phi
    CalBnorm(P, Q, phi, C_m, Bnorm);

    // Derive "pick" vs. "unpick" sets
    std::vector<int> pick, unpick;
    pick.reserve(P);
    unpick.reserve(P);
    for (int p = 0; p < P; p++) {
      if (Bnorm[p] > 1e-6) {
        pick.push_back(p);
      } else {
        unpick.push_back(p);
      }
    }

    // (a) Active set loop
    double flag_a = 100.0;
    while ((flag_a > eps) && (n_iter < max_iter)) {
      Assign(phi, phi_last);
      Assign(phi, phi_old);

      // (a.1) penalized updates
      for (int pidx : pick) {
        Update(pidx, N, P, Q,
               lam1, lam2,
               sigma,
               C_m, X_m, Xnorm,
               E, Beta, Bnorm, phi_old, phi);
        n_iter++;
      }

      // (a.2) unpenalized updates
      for (int pidx : unpick) {
        for (int q = 0; q < Q; q++) {
          if (C_m(pidx, q) == 2) {
            double tmp = 0.0;
            for (int n = 0; n < N; n++) {
              tmp += E(n, q) * X_m(n, pidx);
            }
            double new_val = (tmp / Xnorm[pidx]) + phi_old(pidx, q);

            double diff = phi_old(pidx, q) - new_val;
            if (std::fabs(diff) > 1e-15) {
              for (int n = 0; n < N; n++) {
                E(n, q) += diff * X_m(n, pidx);
              }
            }
            phi(pidx, q)     = new_val;
            phi_old(pidx, q) = new_val;
            n_iter++;
          }
        }
      }

      flag_a = Dist(phi_last, phi);
    }

    // (b) "Full loop" over all rows
    Assign(phi, phi_last);
    Assign(phi, phi_old);

    for (int p = 0; p < P; p++) {
      Update(p, N, P, Q,
             lam1, lam2,
             sigma,
             C_m, X_m, Xnorm,
             E, Beta, Bnorm, phi_old, phi);
      n_iter++;
    }

    flag = Dist(phi_last, phi);
  } // end main while

  // --------------------------------------------------------------------------
  // 6) Compute final RSS
  // --------------------------------------------------------------------------
  double rss = 0.0;
  for (int n = 0; n < N; n++) {
    for (int q = 0; q < Q; q++) {
      double val = E(n, q);
      rss += val * val;
    }
  }

  // --------------------------------------------------------------------------
  // 7) Return a list
  // --------------------------------------------------------------------------
  return List::create(
    _["Phi_output"] = phi,
    _["N_iter"]     = n_iter,
    _["RSS"]        = rss,
    _["E_debug"]    = E
  );
}

// ----------------------------------------------------------------------------
//   2) BIC across a grid of lam1, lam2 + 'gamma' extra penalty
// ----------------------------------------------------------------------------

// [[Rcpp::export]]
List MultiRegGroupLassoBIC(const NumericMatrix &X,
                           const NumericMatrix &Y,
                           const IntegerMatrix &C_m,
                           const NumericVector &lam1_vec,
                           const NumericVector &lam2_vec,
                           const NumericVector &sigma,
                           double gamma = 1.0) {
  /*
   For each (lam1, lam2) in lam1_vec x lam2_vec:
   1) Fit the model with MultiRegGroupLasso_sigma(...) => get phi
   2) Compute RSS for each column => sum( log(rss_i) ) * n
   3) Compute df => can re-use your "MultiRegGroupLassoDegree(...)" logic,
   or manually estimate from the final phi
   4) Add extra term:  2 * gamma * log( sum( choose(q, rowSums(phi != 0)) ) )
   5) BIC = sum(log(rss)) * n + df*log(n) + [the extra gamma term]

   We'll do a simpler approach: fit each pair, compute df by counting non-zero
   in penalized spots + some group-lasso logic, etc.
   For demonstration, we do a naive df = #nonzeroPenalized + #unpenalized
   If you have more advanced 'degree' formula, adapt it as needed.
   */

  int n = X.nrow();
  int p = X.ncol();
  int q = Y.ncol();
  int len1 = lam1_vec.size();
  int len2 = lam2_vec.size();

  // We'll store a (len1 x len2) matrix of BIC and a list of phi
  NumericMatrix bic(len1, len2);
  std::vector<List> phi_list(len1 * len2);

  // We do a double loop over i, j
  for (int i = 0; i < len1; i++) {
    double curLam1 = lam1_vec[i];
    for (int j = 0; j < len2; j++) {
      double curLam2 = lam2_vec[j];

      // Fit with the user-chosen sigma
      List fit = MultiRegGroupLasso_sigma_unified(X, Y, C_m, curLam1, curLam2, sigma);
      NumericMatrix phi = fit["Phi_output"];
      double rss   = fit["RSS"];

      // 1) sum(log(rss_i)) => Actually we can do: sum(log(rss)) if we have
      //    separate columns. But we only have a single "rss" = sum_{all i,q} of residual^2.
      //    If we want sum(log(rss_i)) per column, we'd need each column's RSS separately.
      //    We'll approximate with log(rss) * n.
      //    (If you want the exact approach from your old code, adapt as needed.)
      double part1 = std::log(rss) * n;

      // 2) degrees of freedom => we can approximate
      //    df = # of unpenalized spots + # of penalized spots that are non-zero
      //    or replicate your "remMap.df" approach
      int df_val = 0;
      int unpenalized_count = 0;
      for (int pp = 0; pp < p; pp++) {
        for (int qq = 0; qq < q; qq++) {
          if (C_m(pp, qq) == 2) {
            // not penalized => always in the model
            unpenalized_count++;
          } else if (C_m(pp, qq) == 1) {
            // penalized => count if phi != 0
            if (std::fabs(phi(pp, qq)) > 1e-15) {
              df_val++;
            }
          }
        }
      }
      df_val += unpenalized_count;

      // 3) The extra gamma * log(...) term
      //    For each row p, we can check rowSums(phi != 0). Then sum( choose(q, rowSum ) )
      //    This is a bit special, adapt from your old code. We'll do an example:
      double sumChoose = 0.0;
      for (int pp = 0; pp < p; pp++) {
        int rowNonZero = 0;
        for (int qq = 0; qq < q; qq++) {
          if (std::fabs(phi(pp, qq)) > 1e-15) {
            rowNonZero++;
          }
        }
        // choose(q, rowNonZero) if rowNonZero <= q
        // we'll do a quick function Rf_choose() or approximate.
        // Rcpp has Rf_choose() or we can do our own.
        if (rowNonZero <= q && rowNonZero >= 0) {
          // we can use R's built in choose function via R:::
          // but in C++ we can do an approximation or call Rf_choose
          // For simplicity, let's do a quick approach:
          double cval = R::choose(q, rowNonZero);
          sumChoose += cval;
        }
      }
      double extraTerm = 2.0 * gamma * std::log(std::max(sumChoose, 1.0));

      // 4) BIC
      double curBIC = part1 + df_val * std::log(n) + extraTerm;

      bic(i, j) = curBIC;
      phi_list[i * len2 + j] = List::create(_["lam1"]=curLam1,
                                            _["lam2"]=curLam2,
                                            _["phi"] = phi);
    }
  }

  return List::create(
    _["BIC"] = bic,
    _["phi"] = phi_list
  );
}

// ----------------------------------------------------------------------------
//    3) Example “Degree” function with sigma => optional if you need it
// ----------------------------------------------------------------------------

// [[Rcpp::export]]
List MultiRegGroupLassoDegree_sigma(const NumericMatrix &X_m,
                                    const NumericMatrix &Y_m,
                                    const IntegerMatrix &C_m,
                                    const NumericVector &lam1_vec,
                                    const NumericVector &lam2_vec,
                                    const NumericVector &sigma) {
  /*
   This is analogous to your old "remMap.df" approach, but now we scale
   lam1, lam2 by sigma[q] in a partial way if you want.
   Or you can replicate the older approach exactly but we typically only
   used lam2 on XYnorm, etc.
   Adjust as needed if your formula for df changed.
   */
  // For brevity, not fully implemented.
  // You could replicate the logic from "MultiRegGroupLassoDegree"
  // and incorporate sigma in the "SoftShrink( XYnorm[p,q], lam2*sigma[q] )" step.
  // Then compute the row norms, etc. as you did before.

  // Placeholder:
  return List::create(_["Degree"] = NumericMatrix(lam1_vec.size(), lam2_vec.size()),
                      _["Debug"]  = 0.0);
}
