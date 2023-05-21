#pragma once
#include <cstdlib>
const double G = 6.67e-11;
const double PI = 3.1415926535;

inline double uniform() { return (double)rand() / RAND_MAX; }
inline size_t discrete_uniform(size_t n) { return rand() % n; }

namespace config {
const double dt = 1e-2;
const double t_end = 100;
const double draw_dt = 1;
const size_t canvas_width = 400, canvas_height = 400;
} // namespace config
