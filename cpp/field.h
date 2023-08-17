#pragma once

#include <godot_cpp/classes/mesh_instance3d.hpp>
#include <Eigen/CXX11/Tensor>

namespace godot {

class Field : public MeshInstance3D {
    GDCLASS(Field, MeshInstance3D)
public:

    int n;
    Eigen::Tensor<double, 3> xt, yt, zt, xy, yz, zx;
   
    static void _bind_methods();

    Field();

    void _ready();
    void _process(double delta);
};

} // namespace godot
