#include "field.h"
#include <godot_cpp/classes/array_mesh.hpp>
#include <godot_cpp/core/class_db.hpp>

using namespace godot;
using namespace Eigen;

void Field::_bind_methods() {
  ClassDB::bind_method(D_METHOD("get_N"), &Field::get_N);
  ClassDB::bind_method(D_METHOD("set_N", "p_N"), &Field::set_N);
  ClassDB::add_property("Field", PropertyInfo(Variant::INT, "N"), "set_N",
                        "get_N");
  ClassDB::bind_method(D_METHOD("get_nu"), &Field::get_nu);
  ClassDB::bind_method(D_METHOD("set_nu", "p_nu"), &Field::set_nu);
  ClassDB::add_property(
      "Field",
      PropertyInfo(Variant::FLOAT, "nu", PROPERTY_HINT_RANGE, "0,1,0.01"),
      "set_nu", "get_nu");
}

Field::Field() {
  std::cout << "start constructor" << std::endl;

  Ref<ArrayMesh> mesh;
  mesh.instantiate();
  set_mesh(mesh);

  stokes.b[0](stokes.N / 2, stokes.N / 2, stokes.N / 2) = 36;
  u::print("finish constructor");
}

void Field::_ready() {
  stokes.solve_fft();
  draw();
}

void Field::draw() {
  std::vector arrow = {
      Vector3(1. / 4, 0, 0),        Vector3(3. / 4, 0, 0),
      Vector3(2. / 3, 1. / 24, 0),  Vector3(3. / 4, 0, 0),
      Vector3(2. / 3, -1. / 24, 0), Vector3(3. / 4, 0, 0),
  };
  std::vector vortex = {
      Vector3(1. / 5, 2. / 5, 0), Vector3(3. / 5, 2. / 5, 0),
      Vector3(3. / 5, 1. / 5, 0), Vector3(3. / 5, 3. / 5, 0),
      Vector3(4. / 5, 3. / 5, 0), Vector3(2. / 5, 3. / 5, 0),
      Vector3(2. / 5, 4. / 5, 0), Vector3(2. / 5, 2. / 5, 0),
  };
  std::vector plus = {
      Vector3(2. / 5, 1. / 2, 1. / 2), Vector3(3. / 5, 1. / 2, 1. / 2),
      Vector3(1. / 2, 2. / 5, 1. / 2), Vector3(1. / 2, 3. / 5, 1. / 2),
      Vector3(1. / 2, 1. / 2, 2. / 5), Vector3(1. / 2, 1. / 2, 3. / 5),
  };

  Vector3 axis = Vector3(1, 1, 1).normalized();
  Transform3D flip;
  flip.basis[0][0] = -1;
  flip.origin[0] = 1;

  Array arrays;
  arrays.resize(Mesh::ARRAY_MAX);
  PackedColorArray colors;
  PackedVector3Array verts;

  for (int x = 0; x < stokes.N; x++) {
    for (int y = 0; y < stokes.N; y++) {
      for (int z = 0; z < stokes.N; z++) {
        auto draw_shape = [&](auto &field, int axis_idx,
                              std::vector<Vector3> shape, double color) {
          double str = field(x, y, z);
          Transform3D trans(Basis(axis, Math_TAU * axis_idx / 3),
                            Vector3(x, y, z));
          if (str < 0) {
            trans *= flip;
          }
          str = str * str * color;
          Color c(str, abs(str) - 0.5, -str);

          for (Vector3 vert : shape) {
            colors.append(c);
            verts.append(trans.xform(vert));
          }
        };

        for (int i : axes) {
          draw_shape(stokes.u[i], i + 1, vortex, 1);
          draw_shape(stokes.p, i, plus, -1);
        }
      }
    }
  }

  arrays[Mesh::ARRAY_COLOR] = colors;
  arrays[Mesh::ARRAY_VERTEX] = verts;
  ArrayMesh *mesh = Object::cast_to<ArrayMesh>(*get_mesh());
  mesh->clear_surfaces();
  mesh->add_surface_from_arrays(Mesh::PRIMITIVE_LINES, arrays);
}
