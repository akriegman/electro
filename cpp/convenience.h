#pragma once

#include <Eigen/CXX11/Tensor>
#include <godot_cpp/variant/utility_functions.hpp>

using u = godot::UtilityFunctions;

namespace Eigen {

using Spectrum = Tensor<std::complex<double>, 3>;
using Scalar = Tensor<double, 3>;
using Triple = DSizes<int, 3>;

Spectrum fft(Scalar f) {
  return f.fft<BothParts, FFT_FORWARD>(Triple{0, 1, 2});
}

Scalar ifft(Spectrum f) {
  return f.fft<RealPart, FFT_REVERSE>(Triple{0, 1, 2});
}

template <typename T, int N>
auto antirotate(const Tensor<T, N> &data, const Index dist, const Index dim) {
  auto singleton = [&](Index idx) {
    DSizes<Index, N> out;
    out.fill(0);
    out[dim] = idx;
    return out;
  };
  auto deprextent = [&](Index ext) {
    DSizes<Index, N> out = data.dimensions();
    out[dim] = ext;
    return out;
  };
  if (dist > 0) {
    return (-1 *
            data.slice(singleton(data.dimension(dim) - dist), deprextent(dist)))
        .concatenate(1 * data.slice(singleton(0),
                                    deprextent(data.dimension(dim) - dist)),
                     dim);
  } else {
    return (1 * data.slice(singleton(-dist),
                           deprextent(data.dimension(dim) + dist)))
        .concatenate(-1 * data.slice(singleton(0), deprextent(-dist)), dim);
  }
}

// TODO: find a way to do more antiperiodic operations elegantly
// maybe a signed_ref type which implements +=, *, =, etc. with appropriate
// parity. would also be nice to make this generic in N
template <typename T>
void peq_ap(Tensor<T, 3> &data, T val, const godot::Vector3i &idxs) {
  int parity = 0;
  Vector3i reduced;
  for (int i : axes) {
    parity += idxs[i] / data.dimension(i);
    reduced[i] = u::posmod(idxs[i], data.dimension(i));
  }
  data(reduced[0], reduced[1], reduced[2]) += val * (parity % 2 ? -1 : 1);
}

} // namespace Eigen
