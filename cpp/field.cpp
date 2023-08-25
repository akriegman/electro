#include "field.h"
#include "convenience.h"
#include <godot_cpp/classes/array_mesh.hpp>
#include <godot_cpp/core/class_db.hpp>

using namespace godot;
using namespace Eigen;

void Field::_bind_methods() {}

Field::Field() : n(8), c(0.5), dt(1. / 2), xt(n, n, n) {
  xt.setConstant(0);
  yt = xt.constant(0);
  zt = xt.constant(0);
  xy = xt.constant(0);
  yz = xt.constant(0);
  zx = xt.constant(0);
  xt(4, 4, 4) = 10;
}

void Field::_ready() { time_since_tick = 0; }

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

  Vector3 axis = Vector3(1, 1, 1).normalized();
  Transform3D flip;
  flip.basis[0][0] = -1;
  flip.origin[0] = 1;

  Array arrays;
  arrays.resize(Mesh::ARRAY_MAX);
  PackedColorArray colors;
  PackedVector3Array verts;

  for (int x = 0; x < n; x++) {
    for (int y = 0; y < n; y++) {
      for (int z = 0; z < n; z++) {
        auto draw_shape = [&](auto &field, int axis_idx,
                              std::vector<Vector3> shape, int color) {
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

        draw_shape(xt, 0, arrow, 1);
        draw_shape(yt, 1, arrow, 1);
        draw_shape(zt, 2, arrow, 1);
        draw_shape(xy, 0, vortex, -1);
        draw_shape(yz, 1, vortex, -1);
        draw_shape(zx, 2, vortex, -1);
      }
    }
  }

  arrays[Mesh::ARRAY_COLOR] = colors;
  arrays[Mesh::ARRAY_VERTEX] = verts;
  Ref<ArrayMesh> mesh;
  mesh.instantiate();
  mesh->add_surface_from_arrays(Mesh::PRIMITIVE_LINES, arrays);
  set_mesh(mesh);
}

void Field::_process(double delta) {
  time_since_tick += delta;

  while (time_since_tick >= dt) {
    time_since_tick -= dt;

    // curl E + d/dt B = 0
    xt += (zx - xy - antirotate(zx, 1, 2) + antirotate(xy, 1, 1)) * c * c;
    yt += (xy - yz - antirotate(xy, 1, 0) + antirotate(yz, 1, 2)) * c * c;
    zt += (yz - zx - antirotate(yz, 1, 1) + antirotate(zx, 1, 0)) * c * c;
    // curl H - d/dt D = J
    xy += xt - yt - antirotate(xt, -1, 1) + antirotate(yt, -1, 0);
    yz += yt - zt - antirotate(yt, -1, 2) + antirotate(zt, -1, 1);
    zx += zt - xt - antirotate(zt, -1, 0) + antirotate(xt, -1, 2);
    // note that aliasing is not an issue because the first three
    // are actually happening on the open tick

    draw();
  }
}
