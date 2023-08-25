#pragma once

#include <Eigen/CXX11/Tensor>
#include <initializer_list>

namespace Eigen {

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

// antiperiodic getter
// unfinished, unneeded
template <typename T, int N>
T get_ap(Tensor<T, N> data, std::initializer_list<Index> idxs) {
  int parity = 0;
  DSizes<Index, N> reduced;
  for (int i = 0; i < N; i++) {
    parity += idxs.begin()[i] / data.dimension(i);
    reduced[i] = idxs.begin()[i] % data.dimension(i);
  }
  return data(reduced) * (parity % 2 == 0 ? 1 : -1);
}

} // namespace Eigen
