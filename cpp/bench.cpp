// parts copied from
// https://stackoverflow.com/questions/22387586/measuring-execution-time-of-a-function-in-c

#include <chrono>
#include <iostream>
#include <utility>

#include "stokes.h"

using namespace std;
using namespace Eigen;

typedef chrono::high_resolution_clock::time_point TimeVar;

#define duration(a) chrono::duration_cast<std::chrono::nanoseconds>(a).count()
#define timeNow() chrono::high_resolution_clock::now()

template <typename F, typename... Args> double bench(F func, Args &&...args) {
  TimeVar t1 = timeNow();
  func(forward<Args>(args)...);
  return duration(timeNow() - t1);
}

int main() {
  for (int N = 4; N <= 16; N++) {
    Stokes s(N);
    cout << "==== N = " << N << " ====" << endl;
    cout << "setup:	" << endl;
    cout << "  fft: " << bench([&]() { s.setup_fft(); }) << endl;
    cout << "   qr: " << bench([&]() { s.setup_qr(); }) << endl;
    cout << "solve: " << endl;
    cout << "  fft: " << bench([&]() { s.solve_fft(); }) << endl;
    cout << "   qr: " << bench([&]() { s.solve_qr(); }) << endl;
  }
}
