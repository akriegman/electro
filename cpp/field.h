#pragma once

#include "stokes.h"
#include <godot_cpp/classes/mesh_instance3d.hpp>
#include <utility>

namespace godot {

class Field : public MeshInstance3D {
  GDCLASS(Field, MeshInstance3D)
public:
  Eigen::Stokes stokes;

  static void _bind_methods();

  Field();

  void draw();
  void _ready() override;

  void set_N(const int N) { stokes = Eigen::Stokes(N); }
  int get_N() const { return stokes.N; }
  void set_nu(const double nu) { stokes = Eigen::Stokes(stokes.N, nu); }
  double get_nu() const { return stokes.nu; }
};

} // namespace godot
