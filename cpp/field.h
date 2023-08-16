#pragma once

#include <godot_cpp/classes/mesh_instance3d.hpp>

namespace godot {

class Field : public MeshInstance3D {
    GDCLASS(Field, MeshInstance3D)
public:
   
    static void _bind_methods();

    Field();

    void _process(double delta);
};

}
