#pragma once

#include "convenience.h"
#include <Eigen/Sparse>

namespace Eigen {

class Stokes {
public:
  // ==== parameters ====
  int N;
  double nu;

  // ==== data ====
  // p is the negative pressure, because that convention is better and
  // makes things symmetric.
  Tensor<double, 3> b[3], u[3], p;

  // ==== setup ====
  // U[i] = S[i][j] B[j]
  // P    = S[3][j] B[j]
  Spectrum S[4][3];
  Spectrum denom;
  std::unique_ptr<SparseMatrix<double>> stokes;
  std::unique_ptr<SparseQR<SparseMatrix<double>, AMDOrdering<int>>> solver;

  // ==== methods ====
  Stokes(int N = 8, double nu = 1);

  void setup_fft();
  void setup_qr();
  void solve_fft();
  void solve_qr();
};

} // namespace Eigen
