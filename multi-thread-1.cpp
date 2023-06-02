#include <thread>

#include "algorithms.hpp"

static void multi_thread_1_aux(Scenario &bodies, std::vector<Vect> &curr,
                               std::vector<Vect> &next, size_t start,
                               size_t end) {
    size_t n = bodies.r.size();
    for (size_t i = start; i < end; i++) {
        for (size_t j = i + 1; j < end; ++j) {
            Vect dr = curr[j] - curr[i];

            double dist_sq = std::max(dr.norm2(), 1e-6);
            double dist = sqrt(dist_sq);
            double force_mag = G * bodies.m[i] * bodies.m[j] / dist_sq;
            Vect force = dr * (force_mag / dist);

            bodies.v[i] += force * (dt / bodies.m[i]);
            bodies.v[j] -= force * (dt / bodies.m[j]);
        }

        for (size_t j = 0; j < start; ++j) {
            Vect dr = curr[j] - curr[i];

            double dist_sq = std::max(dr.norm2(), 1e-6);
            double dist = sqrt(dist_sq);
            double force_mag = G * bodies.m[i] * bodies.m[j] / dist_sq;
            Vect force = dr * (force_mag / dist);

            bodies.v[i] += force * (dt / bodies.m[i]);
        }

        for (size_t j = end; j < n; ++j) {
            Vect dr = curr[j] - curr[i];

            double dist_sq = std::max(dr.norm2(), 1e-6);
            double dist = sqrt(dist_sq);
            double force_mag = G * bodies.m[i] * bodies.m[j] / dist_sq;
            Vect force = dr * (force_mag / dist);

            bodies.v[i] += force * (dt / bodies.m[i]);
        }

        next[i] = curr[i] + bodies.v[i] * dt;
    }
}

#ifdef VISUALIZE
void multi_thread_1(Scenario &bodies, size_t n_threads, Drawer &drawer) {
#else
void multi_thread_1(Scenario &bodies, size_t n_threads) {
#endif

    std::cout << "Simulation with multi_thread_1\n";

    size_t n = bodies.r.size();
    size_t chunk_size = n / n_threads;
    std::thread threads[n_threads - 1];
    std::vector<Vect> tmp_r(n);
    std::vector<Vect> *curr = &bodies.r, *next = &tmp_r;

    for (double t = 0; t <= t_end; t += dt) {
#ifdef VISUALIZE
        drawer.trigger_draw(t, curr);
#endif
        for (size_t i = 0; i < n_threads - 1; ++i) {
            threads[i] = std::thread(multi_thread_1_aux, std::ref(bodies),
                                     std::ref(*curr), std::ref(*next),
                                     i * chunk_size, (i + 1) * chunk_size);
        }
        multi_thread_1_aux(bodies, *curr, *next, (n_threads - 1) * chunk_size,
                           n);
        for (size_t i = 0; i < n_threads - 1; ++i) {
            threads[i].join();
        }
        std::swap(curr, next);
    }
}
