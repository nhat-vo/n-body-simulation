#include <thread>

#include "algorithms.hpp"

inline static size_t getForceIndex(size_t body_id, size_t thread_id,
                                   size_t n_threads, size_t n_bodies) {
    // return body_id * n_threads + thread_id; // cache unfriendly
    return body_id + thread_id * n_bodies;
}

static void aux(const Scenario &bodies, std::vector<Vect> &delta_v, size_t load,
                size_t thread_id, size_t n_threads) {
    size_t n = bodies.r.size();
    const auto &r = bodies.r;
    const auto &m = bodies.m;

    size_t mid = n / 2;
    size_t start = thread_id * load;
    size_t end = std::min(start + load, n);

    // assume delta_v is organized as contiguous blocks for each thread
    std::fill(delta_v.begin() + getForceIndex(0, thread_id, n_threads, n),
              delta_v.begin() + getForceIndex(n, thread_id, n_threads, n),
              Vect(0, 0));

    for (size_t i = start; i < end; i++) {
        for (size_t j = std::max(i + 1, mid); j < n; ++j) {
            Vect dr = r[j] - r[i];
            double dist_sq = std::max(dr.norm2(), 1e-6);
            double dist = sqrt(dist_sq);
            double force_mag = G * m[i] * m[j] / dist_sq;
            Vect accel = dr * (force_mag / dist * dt);

            // std::cout << "i: " << i << " j: " << j << std::endl;
            delta_v[getForceIndex(i, thread_id, n_threads, n)] += accel / m[i];
            delta_v[getForceIndex(j, thread_id, n_threads, n)] -= accel / m[j];
        }
    }

    start = (n_threads - 1 - thread_id) * load;
    end = std::min(start + load, n);
    for (size_t i = start; i < end; i++) {
        for (size_t j = i + 1; j < mid; j++) {
            Vect dr = r[j] - r[i];
            double dist_sq = std::max(dr.norm2(), 1e-6);
            double dist = sqrt(dist_sq);
            double force_mag = G * m[i] * m[j] / dist_sq;
            Vect accel = dr * (force_mag / dist * dt);

            // std::cout << "i: " << i << " j: " << j << std::endl;
            delta_v[getForceIndex(i, thread_id, n_threads, n)] += accel / m[i];
            delta_v[getForceIndex(j, thread_id, n_threads, n)] -= accel / m[j];
        }
    }
}

#ifdef VISUALIZE
void multi_thread_2(Scenario &bodies, size_t n_threads, Drawer &drawer) {
#elif defined(WRITE)
void multi_thread_2(Scenario &bodies, size_t n_threads, Writer &writer) {
#else
void multi_thread_2(Scenario &bodies, size_t n_threads) {
#endif

#ifdef DEBUG
    std::cout << "Simulation with multi_thread_2\n";
#endif

    size_t n = bodies.r.size();
    size_t chunk_size = (n + n_threads - 1) / n_threads;
    std::thread threads[n_threads - 1];

    std::vector<Vect> delta_v(n * n_threads);

    for (double t = 0; t < t_end; t += dt) {
#ifdef VISUALIZE
        drawer.trigger_draw(t, &bodies.r);
#elif WRITE
        writer.write(bodies.r);
#endif

        // std::cout << "update\n";
        // std::fill(delta_v.begin(), delta_v.end(), Vect(0, 0));
        for (size_t i = 0; i < n_threads - 1; ++i) {
            threads[i] = std::thread(aux, std::ref(bodies), std::ref(delta_v),
                                     chunk_size, i, n_threads);
        }
        aux(bodies, delta_v, chunk_size, n_threads - 1, n_threads);
        for (size_t i = 0; i < n_threads - 1; ++i) {
            threads[i].join();
        }

        for (size_t j = 0; j < n_threads; j++) {
            for (size_t i = 0; i < n; i++) {
                bodies.v[i] += delta_v[getForceIndex(i, j, n_threads, n)];
            }
        }
        for (size_t i = 0; i < n; i++) bodies.r[i] += bodies.v[i] * dt;
    }
}
