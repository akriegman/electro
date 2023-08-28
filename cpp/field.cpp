#include "field.h"
#include "convenience.h"
#include <godot_cpp/classes/array_mesh.hpp>
#include <godot_cpp/core/class_db.hpp>

using namespace godot;
using namespace Eigen;

void Field::_bind_methods() {
  ClassDB::bind_method(D_METHOD("get_N"), &Field::get_N);
  ClassDB::bind_method(D_METHOD("set_N", "p_N"), &Field::set_N);
  ClassDB::add_property("Field", PropertyInfo(Variant::INT, "N"), "set_N",
                        "get_N");
  ClassDB::bind_method(D_METHOD("get_c"), &Field::get_c);
  ClassDB::bind_method(D_METHOD("set_c", "p_c"), &Field::set_c);
  ClassDB::add_property(
      "Field",
      PropertyInfo(Variant::FLOAT, "c", PROPERTY_HINT_RANGE, "0,1,0.01"),
      "set_c", "get_c");
  ClassDB::bind_method(D_METHOD("get_dt"), &Field::get_dt);
  ClassDB::bind_method(D_METHOD("set_dt", "p_dt"), &Field::set_dt);
  ClassDB::add_property(
      "Field",
      PropertyInfo(Variant::FLOAT, "dt", PROPERTY_HINT_RANGE, "0,1,0.01"),
      "set_dt", "get_dt");
  ClassDB::bind_method(D_METHOD("make_charge", "pos"), &Field::make_charge);
  ClassDB::bind_method(D_METHOD("charge_moved", "pos"), &Field::charge_moved);
}

int n_threads() { return std::thread::hardware_concurrency() - 1; }

Field::Field()
    : N(8), c(0.5), dt(1. / 2), q(6), tp(n_threads()), tpd(&tp, n_threads()) {

  Tensor<double, 3> empty(N, N, N);
  empty.setZero();
  for (int i : axes) {
    H[i] = D[i] = J[i] = empty;
  }

  Ref<ArrayMesh> mesh;
  mesh.instantiate();
  set_mesh(mesh);
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

  for (int x = 0; x < N; x++) {
    for (int y = 0; y < N; y++) {
      for (int z = 0; z < N; z++) {
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
          draw_shape(H[i], i, arrow, 1);
          draw_shape(D[i], i, vortex, -c);
          // draw_shape(J[i], (i + 1) % 3, vortex, 1);
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

void Field::_process(double delta) {
  time_since_tick += delta;

  while (time_since_tick >= dt) {
    time_since_tick -= dt;

    // closed tick
    // curl H - d/dt D = J
    for (int i : axes) {
      int j = (i + 1) % 3;
      int k = (i + 2) % 3;
      D[i].device(tpd) += H[i] - H[j] - antirotate(H[i], -1, j) +
                          antirotate(H[j], -1, i) - J[k];
      J[k].setZero();
    }
    // open tick
    // curl E + d/dt B = 0
    for (int i : axes) {
      int j = (i + 1) % 3;
      int k = (i + 2) % 3;
      H[i].device(tpd) +=
          (D[k] - D[i] - antirotate(D[k], 1, k) + antirotate(D[i], 1, j)) * c *
          c;
    }
    // note that aliasing is not an issue because of how we split up the updates

    draw();
  }
}

void Field::make_charge(Vector3 pos) {
  charge_cell = pos.floor();
  for (int i : axes) {
    Vector3i face = charge_cell;
    face[i]++;
    peq_ap(D[i], -q / 6, charge_cell);
    peq_ap(D[i], q / 6, face);
  }
}

void Field::charge_moved(Vector3 pos) {
  Vector3i new_cell = pos.floor();

  for (int i : axes) {
    // at most one of these loops runs
    while (charge_cell[i] < new_cell[i]) {
      charge_cell[i]++;
      peq_ap(J[i], q, charge_cell);
    }
    while (charge_cell[i] > new_cell[i]) {
      peq_ap(J[i], -q, charge_cell);
      charge_cell[i]--;
    }
  }
}
