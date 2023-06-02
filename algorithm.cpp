// #include "common.hpp"
// #include "visualizer.hpp"
#include <chrono>
#include <cmath>
#include <iostream>
#include <list>
#include <thread>
#include <vector>

#include "barnes-hut.hpp"

using namespace config;

void compute_forces(Scenario &bodies, size_t start, size_t end, double dt) {
    size_t n = bodies.r.size();
    for (size_t i = start; i < end; i++) {
        for (size_t j = i + 1; j < n; j++) {
            Vect dr = bodies.r[j] - bodies.r[i];

            double dist_sq = std::max(dr.norm2(), 1e-6);
            double dist = sqrt(dist_sq);
            double force_mag = G * bodies.m[i] * bodies.m[j] / dist_sq;
            Vect force = dr * (force_mag / dist);

            bodies.v[i] += force * (dt / bodies.m[i]);
            bodies.v[j] -= force * (dt / bodies.m[j]);
        }
    }
}

void update_positions(Scenario &bodies, int start, int end, double dt) {
    for (int i = start; i < end; i++)
        bodies.r[i] = bodies.r[i] + bodies.v[i] * dt;
}

#ifdef VISUALIZE
void single_thread(Scenario &bodies, size_t n_threads, Drawer &drawer) {
#else
void single_thread(Scenario &bodies, size_t n_threads) {
#endif

    std::cout << "Simulation with single_thread\n";

    for (double t = 0; t <= t_end; t += dt) {
#ifdef VISUALIZE
        drawer.trigger_draw(t, &bodies.r);
#endif

        compute_forces(bodies, 0, bodies.r.size(), dt);
        update_positions(bodies, 0, bodies.r.size(), dt);
    }
}

#ifdef VISUALIZE
void barnes_hut(Scenario &bodies, size_t n_threads, Drawer &drawer) {
#else
void barnes_hut(Scenario &bodies, size_t n_threads) {
#endif

    std::cout << "Simulation with single-threaded Barnes-Hut\n";

    for (double t = 0; t <= t_end; t += dt) {
#ifdef VISUALIZE
        drawer.trigger_draw(t, &bodies.r);
#endif

        barnes_hut_update_step(bodies);
    }
}

void multi_thread_1_aux(Scenario &bodies, std::vector<Vect> &curr,
                        std::vector<Vect> &next, size_t start, size_t end) {
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

int main(int argc, char **argv) {
    if (argc == 1) {
        std::cout << "Usage: " << argv[0]
                  << "number of threads> <number of bodies (optional, default "
                     "to 100)>"
                  << std::endl;
        return 0;
    }

    size_t n_threads = 0, n_bodies = 0;
    if (argc >= 2) {
        try {
            int tmp = std::stoi(argv[1]);
            if (tmp <= 0) {
                std::cout << "Number of thread should be at least 1."
                          << std::endl;
                return 0;
            }
            n_threads = tmp;
        } catch (...) {
            std::cout << "Invalid argument for <number of threads>"
                      << std::endl;
        }
    }
    if (argc >= 3) {
        try {
            int tmp = std::stoi(argv[2]);
            if (tmp <= 0) {
                std::cout << "Number of bodies should be at least 1."
                          << std::endl;
                return 0;
            }
            n_bodies = tmp;
        } catch (...) {
            std::cout << "Invalid argument for <number of bodies>" << std::endl;
        }
    } else {
        n_bodies = 100;
    }

#ifdef VISUALIZE
    Magick::InitializeMagick(*argv);
#endif

    // setup scenario
    const Vect offset = {canvas_width / 2, canvas_height / 2};
    // double v = sqrt(1 * G * 1e14 / 50);
    // Scenario bodies{
    //     {1e14, 1},
    //     {offset + Vect(0, 0), offset + Vect(0, -50)},
    //     {{0, 0}, {v, 0}},
    //     {"red", "green"},
    // };

    Scenario bodies{{1e14}, {offset + Vect(0, 0)}, {{0, 0}}, {"red"}};
    std::vector<std::string> colors{"blue", "green",  "gold",   "grey",
                                    "pink", "orange", "purple", "brown"};
    for (size_t i = 0; i < n_bodies; ++i) {
        bodies.m.push_back(10 + 5 * uniform());
        bodies.colors.push_back(colors[rand() % colors.size()]);
        bodies.r.push_back(Vect((0.25 + 0.5 * uniform()) * canvas_width,
                                (0.25 + 0.5 * uniform()) * canvas_height));

        Vect dir = bodies.r.back() - bodies.r.front();
        double dist = dir.norm();
        double v = sqrt((0.5 + uniform()) * G * bodies.m.front() / dist);
        dir = dir / dist;
        bodies.v.push_back({v * (-dir.y), v * dir.x});
    }

#ifdef VISUALIZE
    // setup drawing thread
    Drawer drawer(bodies.colors);
#endif

    std::cout << "Starting simulation of " << n_bodies << " bodies with "
              << n_threads << " parallel threads..." << std::endl;

    auto start = std::chrono::steady_clock::now();
#ifdef VISUALIZE
    barnes_hut(bodies, n_threads, drawer);
#else
    // multi_thread_1(bodies, n_threads);
    single_thread(bodies, n_threads);
#endif
    auto end = std::chrono::steady_clock::now();

    std::chrono::duration<double> duration = end - start;
    std::cout << "Simulation ended. Took " << duration.count() << "s"
              << std::endl;
}
