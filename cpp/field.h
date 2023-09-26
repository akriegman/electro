#pragma once

#define EIGEN_USE_THREADS
#include <Eigen/CXX11/Tensor>
#include <godot_cpp/classes/mesh_instance3d.hpp>

const int axes[] = {0, 1, 2};

namespace godot {

class Field : public MeshInstance3D {
  GDCLASS(Field, MeshInstance3D)
public:
  int N;
  // p is the negative pressure, because that convention is better and
  // makes things symmetric.
  Eigen::Tensor<double, 3> b[3], u[3], p;
  double nu;

  static void _bind_methods();

  Field();

  void draw();
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
