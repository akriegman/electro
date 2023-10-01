#include "field.h"
#include "convenience.h"
#include <Eigen/Sparse>
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

Field::Field() : N(8), nu(1) {

  Tensor<double, 3> empty(N, N, N);
  empty.setZero();
  for (int i : axes) {
    b[i] = u[i] = p = empty;
  }

  Ref<ArrayMesh> mesh;
  mesh.instantiate();
  set_mesh(mesh);
}

void Field::_ready() {
  b[0](N / 2, N / 2, N / 2) = 24;

  Spectrum B[3];
  for (int i : axes) {
    B[i] = fft(b[i]);
  }

  // U[i] = M[i][j] B[j]
  // P    = M[3][j] B[j]
  Spectrum M[4][3];
  Scalar tmp(N, N, N);

  tmp.setZero();
  tmp(0, 0, 0) = -4;
  tmp(0, 1, 0) = tmp(0, N - 1, 0) = tmp(0, 0, 1) = tmp(0, 0, N - 1) = 1;
  M[0][0] = fft(tmp);
  tmp = tmp.shuffle(Triple{2, 0, 1}).eval();
  M[1][1] = fft(tmp);
  tmp = tmp.shuffle(Triple{2, 0, 1}).eval();
  M[2][2] = fft(tmp);

  tmp.setZero();
  tmp(0, 0, 0) = tmp(1, N - 1, 0) = 1;
  tmp(1, 0, 0) = tmp(0, N - 1, 0) = -1;
  M[0][1] = fft(tmp.shuffle(Triple{0, 1, 2}));
  M[0][2] = fft(tmp.shuffle(Triple{0, 2, 1}));
  M[1][2] = fft(tmp.shuffle(Triple{2, 0, 1}));
  M[1][0] = fft(tmp.shuffle(Triple{1, 0, 2}));
  M[2][0] = fft(tmp.shuffle(Triple{1, 2, 0}));
  M[2][1] = fft(tmp.shuffle(Triple{2, 1, 0}));

  tmp.setZero();
  tmp(0, 0, 0) = 7;
  tmp(N - 1, 0, 0) = -7;
  tmp(0, 1, 0) = tmp(0, N - 1, 0) = tmp(0, 0, 1) = tmp(0, 0, N - 1) =
      tmp(1, 0, 0) = -1;
  tmp(N - 1, 1, 0) = tmp(N - 1, N - 1, 0) = tmp(N - 1, 0, 1) =
      tmp(N - 1, 0, N - 1) = tmp(N - 2, 0, 0) = 1;
  M[3][0] = nu * fft(tmp);
  M[3][1] = nu * fft(tmp.shuffle(Triple{2, 0, 1}));
  M[3][2] = nu * fft(tmp.shuffle(Triple{1, 2, 0}));

  tmp.setZero();
  tmp(0, 0, 0) = -6;
  tmp(1, 0, 0) = tmp(N - 1, 0, 0) = tmp(0, 1, 0) = tmp(0, N - 1, 0) =
      tmp(0, 0, 1) = tmp(0, 0, N - 1) = 1;
  Spectrum denom = nu * fft(tmp).pow(2);

  for (int i : axes) {
    Spectrum U = (M[i][0] * B[0] + M[i][1] * B[1] + M[i][2] * B[2]) / denom;
    U(0, 0, 0) = 0;
    u[i] = ifft(U);
  }
  Spectrum P = (M[3][0] * B[0] + M[3][1] * B[1] + M[3][2] * B[2]) / denom;
  P(0, 0, 0) = 0;
  p = ifft(P);

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
          draw_shape(u[i], i + 1, vortex, 1);
          draw_shape(p, i, plus, -1);
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
