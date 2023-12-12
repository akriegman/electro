#include "stokes.h"
#include <iostream>

using namespace Eigen;
using namespace std;

Stokes::Stokes(int N, double nu) : N(N), nu(nu) {

  Tensor<double, 3> empty(N, N, N);
  empty.setZero();
  for (int i : axes) {
    b[i] = u[i] = empty;
  }
  p = empty;
}

void Stokes::setup_fft() {

  Scalar tmp(N, N, N);

  tmp.setZero();
  tmp(0, 0, 0) = -4;
  tmp(0, 1, 0) = tmp(0, N - 1, 0) = tmp(0, 0, 1) = tmp(0, 0, N - 1) = 1;
  S[0][0] = fft(tmp);
  tmp = tmp.shuffle(Triple{2, 0, 1}).eval();
  S[1][1] = fft(tmp);
  tmp = tmp.shuffle(Triple{2, 0, 1}).eval();
  S[2][2] = fft(tmp);

  tmp.setZero();
  tmp(0, 0, 0) = tmp(1, N - 1, 0) = 1;
  tmp(1, 0, 0) = tmp(0, N - 1, 0) = -1;
  S[0][1] = fft(tmp.shuffle(Triple{0, 1, 2}));
  S[0][2] = fft(tmp.shuffle(Triple{0, 2, 1}));
  S[1][2] = fft(tmp.shuffle(Triple{2, 0, 1}));
  S[1][0] = fft(tmp.shuffle(Triple{1, 0, 2}));
  S[2][0] = fft(tmp.shuffle(Triple{1, 2, 0}));
  S[2][1] = fft(tmp.shuffle(Triple{2, 1, 0}));

  tmp.setZero();
  tmp(0, 0, 0) = 7;
  tmp(N - 1, 0, 0) = -7;
  tmp(0, 1, 0) = tmp(0, N - 1, 0) = tmp(0, 0, 1) = tmp(0, 0, N - 1) =
      tmp(1, 0, 0) = -1;
  tmp(N - 1, 1, 0) = tmp(N - 1, N - 1, 0) = tmp(N - 1, 0, 1) =
      tmp(N - 1, 0, N - 1) = tmp(N - 2, 0, 0) = 1;
  S[3][0] = nu * fft(tmp);
  S[3][1] = nu * fft(tmp.shuffle(Triple{2, 0, 1}));
  S[3][2] = nu * fft(tmp.shuffle(Triple{1, 2, 0}));

  tmp.setZero();
  tmp(0, 0, 0) = -6;
  tmp(1, 0, 0) = tmp(N - 1, 0, 0) = tmp(0, 1, 0) = tmp(0, N - 1, 0) =
      tmp(0, 0, 1) = tmp(0, 0, N - 1) = 1;
  denom = -nu * fft(tmp).pow(2);
}

void Stokes::setup_qr() {
  int M = N * N * N;

  // written with ChatGPT:

  stokes = make_unique<SparseMatrix<double>>(4 * M, 4 * M);
  vector<Triplet<double>> triplets;

  for (int x = 0; x < N; ++x) {
    for (int y = 0; y < N; ++y) {
      for (int z = 0; z < N; ++z) {

        int idx = x + N * (y + N * z);
        int neighbors[] = {
            posmod(x - 1, N) + N * y + N * N * z,
            posmod(x + 1, N) + N * y + N * N * z,
            x + N * posmod(y - 1, N) + N * N * z,
            x + N * posmod(y + 1, N) + N * N * z,
            x + N * y + N * N * posmod(z - 1, N),
            x + N * y + N * N * posmod(z + 1, N),
        };

        for (int i : axes) {
          int off = i * M;
          int pff = 3 * M;

          // -grad p (here p is the negative pressure)
          triplets.push_back(Triplet<double>(pff + idx, off + idx, -1));
          triplets.push_back(
              Triplet<double>(pff + neighbors[2 * i], off + idx, 1));

          // -nu lap u
          triplets.push_back(Triplet<double>(off + idx, off + idx, 6 * nu));
          for (int neighbor : neighbors) {
            triplets.push_back(Triplet<double>(off + neighbor, off + idx, -nu));
          }

          // div u
          triplets.push_back(Triplet<double>(off + idx, pff + idx, -1));
          triplets.push_back(
              Triplet<double>(off + neighbors[2 * i + 1], pff + idx, 1));
        }
      }
    }
  }

  stokes->setFromTriplets(triplets.begin(), triplets.end());

  solver = make_unique<SparseQR<SparseMatrix<double>, AMDOrdering<int>>>();
  solver->compute(*stokes);

  if (solver->info() != Success) {
    cout << "problem decomposing stokes: " << solver->info() << endl;
  }
}

void Stokes::solve_fft() {

  Spectrum B[3];
  for (int i : axes) {
    B[i] = fft(b[i]);
  }

  for (int i : axes) {
    Spectrum U = (S[i][0] * B[0] + S[i][1] * B[1] + S[i][2] * B[2]) / denom;
    U(0, 0, 0) = 0;
    u[i] = ifft(U);
  }
  Spectrum P = (S[3][0] * B[0] + S[3][1] * B[1] + S[3][2] * B[2]) / denom;
  P(0, 0, 0) = 0;
  p = ifft(P);

  // ==== validate ====
  // VectorXd u_p_flat = pack(u, p);
  // Scalar bob[3];
  // Scalar null;
  // copy(begin(b), end(b), begin(bob));
  // unpack(stokes * u_p_flat, bob, null);
  // cout << bob[0] << endl;
  // cout << null << endl;
}

void Stokes::solve_qr() {
  int M = N * N * N;

  VectorXd b_flat = pack(b, b[0].constant(0));
  VectorXd u_p_flat = solver->solve(b_flat);

  if (solver->info() != Success) {
    cout << "problem solving stokes: " << solver->info() << endl;
  }

  unpack(u_p_flat, u, p);

  // cancel net flow
  for (int i : axes) {
    Tensor<double, 0> mean = u[i].mean();
    u[i] -= u[i].constant(mean());
  }
  // pressure normalization is arbitrary
  Tensor<double, 0> mean = p.mean();
  p -= p.constant(mean());

  // ==== validate ====

  // VectorXd b_flat_back = stokes * u_p_flat;
  // Scalar bx = TensorMap<Scalar>(&b_flat_back[0], N, N, N);
  // cout << bx << endl;
}
