#pragma once

#include <Eigen/CXX11/Tensor>
#include <godot_cpp/classes/mesh_instance3d.hpp>

const int axes[] = {0, 1, 2};

namespace godot {

class Field : public MeshInstance3D {
  GDCLASS(Field, MeshInstance3D)
public:
  int N;
  Eigen::Tensor<double, 3> H[3], D[3], J[3];
  double c;
  double dt;
  double time_since_tick;
  Vector3i charge_cell;
  double q;

  static void _bind_methods();

  Field();

  void draw();
  void _ready();
  void _process(double delta);

  void make_charge(Vector3 pos);
  void charge_moved(Vector3 pos);

  void set_N(const int p_N) {
    N = p_N;
    for (int i : axes) {
      H[i] = D[i] = J[i] = Eigen::Tensor<double, 3>(N, N, N).setZero();
    }
  }
  int get_N() const { return N; }
  void set_c(const double p_c) { c = p_c; }
  double get_c() const { return c; }
  void set_dt(const double p_dt) { dt = p_dt; }
  double get_dt() const { return dt; }
};

} // namespace godot
