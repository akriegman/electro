#pragma once

#define EIGEN_USE_THREADS
#include "convenience.h"
#include <Eigen/CXX11/Tensor>
#include <Eigen/Sparse>
#include <godot_cpp/classes/mesh_instance3d.hpp>

namespace godot {

class Field : public MeshInstance3D {
  GDCLASS(Field, MeshInstance3D)
public:
  // ==== parameters ====
  int N;
  double nu;

  // ==== data ====
  // p is the negative pressure, because that convention is better and
  // makes things symmetric.
  Eigen::Tensor<double, 3> b[3], u[3], p;

  // ==== setup ====
  // U[i] = S[i][j] B[j]
  // P    = S[3][j] B[j]
  Eigen::Spectrum S[4][3];
  Eigen::Spectrum denom;
  Eigen::SparseMatrix<double> stokes;
  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::AMDOrdering<int>> solver;

  static void _bind_methods();

  Field();

  void draw();
  void solve_fft();
  void solve_qr();
  void _ready() override;

  void set_N(const int p_N) {
    N = p_N;
    for (int i : axes) {
      b[i] = u[i] = p = Eigen::Tensor<double, 3>(N, N, N).setZero();
    }
  }
  int get_N() const { return N; }
  void set_nu(const double p_nu) { nu = p_nu; }
  double get_nu() const { return nu; }
};

} // namespace godot
