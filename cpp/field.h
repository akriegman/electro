#pragma once

#include <Eigen/CXX11/Tensor>
#include <godot_cpp/classes/mesh_instance3d.hpp>

namespace godot {

class Field : public MeshInstance3D {
  GDCLASS(Field, MeshInstance3D)
public:
  int n;
  Eigen::Tensor<double, 3> xt, yt, zt, xy, yz, zx;
  double c;
  double dt;
  double time_since_tick;

  static void _bind_methods();

  Field();

  void draw();
  void _ready();
  void _process(double delta);
};

} // namespace godot
