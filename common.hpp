#pragma once
#include <vector>

#include "vect.hpp"

const double G = 6.67e-11;
const double PI = 3.1415926535;
const double theta = 2.;

inline double uniform() { return (double)rand() / RAND_MAX; }
inline size_t discrete_uniform(size_t n) { return rand() % n; }

namespace config {
const double dt = 1e-1;
const double t_end = 100 + 1e-2;
const double draw_dt = 1;
const int canvas_width = 400, canvas_height = 400;

// TODO: this is a bit of a hard code for Barnes-Hut
const int universe_width = 4000, universe_height = 4000;
}  // namespace config

struct Scenario {
    std::vector<double> m;
    std::vector<Vect> r;
    std::vector<Vect> v;
    std::vector<std::string> colors;
};
