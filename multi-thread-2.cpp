#include <thread>

#include "algorithms.hpp"

inline static size_t getForceIndex(size_t body_id, size_t thread_id,
                                   size_t n_threads) {
    return body_id * n_threads + thread_id;
}

static void aux(const Scenario &bodies, std::vector<Vect> &delta_v,
                size_t start, size_t end, size_t thread_id, size_t n_threads) {
    size_t n = bodies.r.size();
    const auto &r = bodies.r;
    const auto &m = bodies.m;

    for (size_t i = start; i < end; i++) {
        for (size_t j = i + 1; j < n; ++j) {
            Vect dr = r[j] - r[i];
            double dist_sq = std::max(dr.norm2(), 1e-6);
            double dist = sqrt(dist_sq);
            double force_mag = G * m[i] * m[j] / dist_sq;
            Vect dv = dr * (force_mag / dist * dt);

            delta_v[getForceIndex(i, thread_id, n_threads)] += dv / m[i];
            delta_v[getForceIndex(j, thread_id, n_threads)] -= dv / m[j];
        }
    }
}

#ifdef VISUALIZE
void multi_thread_2(Scenario &bodies, size_t n_threads, Drawer &drawer) {
#else
void multi_thread_2(Scenario &bodies, size_t n_threads) {
#endif

    std::cout << "Simulation with multi_thread_2\n";

    size_t n = bodies.r.size();
    size_t chunk_size = n / n_threads;
    std::thread threads[n_threads - 1];

    std::vector<Vect> delta_v(n * n_threads);

    for (double t = 0; t <= t_end; t += dt) {
#ifdef VISUALIZE
        drawer.trigger_draw(t, &bodies.r);
#endif

        std::fill(delta_v.begin(), delta_v.end(), Vect(0, 0));
        for (size_t i = 0; i < n_threads - 1; ++i) {
            threads[i] =
                std::thread(aux, std::ref(bodies), std::ref(delta_v),
                            i * chunk_size, (i + 1) * chunk_size, i, n_threads);
        }
        aux(bodies, delta_v, (n_threads - 1) * chunk_size, n, n_threads - 1,
            n_threads);
        for (size_t i = 0; i < n_threads - 1; ++i) {
            threads[i].join();
        }

        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n_threads; j++) {
                bodies.v[i] += delta_v[getForceIndex(i, j, n_threads)];
            }
            bodies.r[i] += bodies.v[i] * dt;
        }
    }
}
