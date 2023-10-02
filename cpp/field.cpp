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
  // ==== setup fft ====

  Scalar tmp(N, N, N);

  tmp.setZero();
  tmp(0, 0, 0) = -4;
  tmp(0, 1, 0) = tmp(0, N - 1, 0) = tmp(0, 0, 1) = tmp(0, 0, N - 1) = 1;
  S[0][0] = fft(tmp);
  tmp = tmp.shuffle(Triple{2, 0, 1}).eval();
  S[1][1] = fft(tmp);
  tmp = tmp.shuffle(Triple{2, 0, 1}).eval();
  S[2][2] = fft(tmp);

  tmp.setZero();
  tmp(0, 0, 0) = tmp(1, N - 1, 0) = 1;
  tmp(1, 0, 0) = tmp(0, N - 1, 0) = -1;
  S[0][1] = fft(tmp.shuffle(Triple{0, 1, 2}));
  S[0][2] = fft(tmp.shuffle(Triple{0, 2, 1}));
  S[1][2] = fft(tmp.shuffle(Triple{2, 0, 1}));
  S[1][0] = fft(tmp.shuffle(Triple{1, 0, 2}));
  S[2][0] = fft(tmp.shuffle(Triple{1, 2, 0}));
  S[2][1] = fft(tmp.shuffle(Triple{2, 1, 0}));

  tmp.setZero();
  tmp(0, 0, 0) = 7;
  tmp(N - 1, 0, 0) = -7;
  tmp(0, 1, 0) = tmp(0, N - 1, 0) = tmp(0, 0, 1) = tmp(0, 0, N - 1) =
      tmp(1, 0, 0) = -1;
  tmp(N - 1, 1, 0) = tmp(N - 1, N - 1, 0) = tmp(N - 1, 0, 1) =
      tmp(N - 1, 0, N - 1) = tmp(N - 2, 0, 0) = 1;
  S[3][0] = nu * fft(tmp);
  S[3][1] = nu * fft(tmp.shuffle(Triple{2, 0, 1}));
  S[3][2] = nu * fft(tmp.shuffle(Triple{1, 2, 0}));

  tmp.setZero();
  tmp(0, 0, 0) = -6;
  tmp(1, 0, 0) = tmp(N - 1, 0, 0) = tmp(0, 1, 0) = tmp(0, N - 1, 0) =
      tmp(0, 0, 1) = tmp(0, 0, N - 1) = 1;
  denom = nu * fft(tmp).pow(2);

  // ==== setup qr ====

  int M = N * N * N;

  // written with ChatGPT:

  stokes = SparseMatrix<double>(4 * M, 4 * M);
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

  solver.compute(stokes);

  if (solver.info() != Success) {
    u::print("problem decomposing stokes: ", solver.info());
  }

  // ==== do other stuff ====

  b[0](N / 2, N / 2, N / 2) = 36;
  solve_qr();
  draw();
}

void Field::solve_fft() {

  Spectrum B[3];
  for (int i : axes) {
    B[i] = fft(b[i]);
  }

  for (int i : axes) {
    Spectrum U = (S[i][0] * B[0] + S[i][1] * B[1] + S[i][2] * B[2]) / denom;
    U(0, 0, 0) = 0;
    u[i] = ifft(U);
  }
  Spectrum P = (S[3][0] * B[0] + S[3][1] * B[1] + S[3][2] * B[2]) / denom;
  P(0, 0, 0) = 0;
  p = ifft(P);

  // ==== validate ====
  // VectorXd u_p_flat = pack(u, p);
  // Scalar bob[3];
  // Scalar null;
  // std::copy(std::begin(b), std::end(b), std::begin(bob));
  // unpack(stokes * u_p_flat, bob, null);
  // std::cout << bob[0] << std::endl;
  // std::cout << null << std::endl;
}

void Field::solve_qr() {
  int M = N * N * N;

  VectorXd b_flat = pack(b, b[0].constant(0));
  VectorXd u_p_flat = solver.solve(b_flat);

  if (solver.info() != Success) {
    u::print("problem solving stokes: ", solver.info());
  }

  unpack(u_p_flat, u, p);

  // cancel net flow
  for (int i : axes) {
    Tensor<double, 0> mean = u[i].mean();
    u[i] -= u[i].constant(mean());
  }
  // pressure normalization is arbitrary
  Tensor<double, 0> mean = p.mean();
  p -= p.constant(mean());

  // ==== validate ====

  // VectorXd b_flat_back = stokes * u_p_flat;
  // Scalar bx = TensorMap<Scalar>(&b_flat_back[0], N, N, N);
  // std::cout << bx << std::endl;
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
