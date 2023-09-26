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
  b[2](N / 2, N / 2, N / 2) = 24;

  int M = N * N * N;
  VectorXd b_flat(4 * M);
  for (int i : axes) {
    b_flat.segment(i * M, M) = Map<VectorXd>(b[i].data(), b[i].size());
  }
  b_flat.segment(3 * M, M) = VectorXd::Zero(M);

  // written with ChatGPT:

  SparseMatrix<double> stokes(4 * M, 4 * M);
  std::vector<Triplet<double>> triplets;

  for (int x = 0; x < N; ++x) {
    for (int y = 0; y < N; ++y) {
      for (int z = 0; z < N; ++z) {

        int idx = x + N * (y + N * z);
        int neighbors[] = {
            (int)u::posmod(x - 1, N) + N * y + N * N * z,
            (int)u::posmod(x + 1, N) + N * y + N * N * z,
            x + N * (int)u::posmod(y - 1, N) + N * N * z,
            x + N * (int)u::posmod(y + 1, N) + N * N * z,
            x + N * y + N * N * (int)u::posmod(z - 1, N),
            x + N * y + N * N * (int)u::posmod(z + 1, N),
        };

        for (int i : axes) {
          int off = i * M;
          int pff = 3 * M;

          // -grad p (here p is the negative pressure)
          triplets.push_back(Triplet<double>(pff + idx, off + idx, -1));
          triplets.push_back(
              Triplet<double>(pff + neighbors[2 * i], off + idx, 1));

          // -nu lap u
          triplets.push_back(Triplet<double>(off + idx, off + idx, 6 * nu));
          for (int neighbor : neighbors) {
            triplets.push_back(Triplet<double>(off + neighbor, off + idx, -nu));
          }

          // div u
          triplets.push_back(Triplet<double>(off + idx, pff + idx, -1));
          triplets.push_back(
              Triplet<double>(off + neighbors[2 * i + 1], pff + idx, 1));
        }
      }
    }
  }

  stokes.setFromTriplets(triplets.begin(), triplets.end());

  SparseQR<SparseMatrix<double>, AMDOrdering<int>> solver;
  solver.compute(stokes);

  if (solver.info() != Success) {
    u::print("problem decomposing stokes: ", solver.info());
  }

  VectorXd u_p_flat = solver.solve(b_flat);

  if (solver.info() != Success) {
    u::print("problem solving stokes: ", solver.info());
  }

  u::print(b_flat.isApprox(stokes * u_p_flat));

  for (int i : axes) {
    u[i] = TensorMap<Tensor<double, 3>>(&u_p_flat[i * M], N, N, N);
  }
  p = TensorMap<Tensor<double, 3>>(&u_p_flat[3 * M], N, N, N);

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
