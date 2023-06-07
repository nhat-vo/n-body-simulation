#include <chrono>

#include "algorithms.hpp"

int main(int argc, char **argv) {
    // nothing fancy here, just parsing command line arguments
    if (argc == 1) {
        std::cout << "Usage: " << argv[0] << "<algorithm id (default = 0)>"
                  << "<number of threads (default = 1)>"
                  << "<number of bodies (default = 100)>" << std::endl;

        std::cout << "Algorithms:" << std::endl;
        std::cout << "0: single-thread" << std::endl;
        std::cout << "1: multi-thread-1" << std::endl;
        std::cout << "2: multi-thread-2" << std::endl;
        std::cout << "3: barnes-hut" << std::endl;

        return 0;
    }

    size_t n_threads = 1, n_bodies = 100, algo = 3;
    if (argc >= 2) {
        try {
            int tmp = std::stoi(argv[1]);
            if (tmp < 0 || tmp > 4) {
                std::cout << "Unknown algorithm" << std::endl;
                return 1;
            }
            algo = tmp;
        } catch (...) {
            std::cout << "Invalid argument for <algorithm id>" << std::endl;
            return 1;
        }
    }

    if (argc >= 3) {
        try {
            int tmp = std::stoi(argv[2]);
            if (tmp <= 0) {
                std::cout << "Number of threads should be positive."
                          << std::endl;
                return 1;
            }
            n_threads = tmp;
        } catch (...) {
            std::cout << "Invalid argument for <number of threads>"
                      << std::endl;
            return 1;
        }
    }

    if (argc >= 4) {
        try {
            int tmp = std::stoi(argv[3]);
            if (tmp <= 0) {
                std::cout << "Number of bodies should be at least 1."
                          << std::endl;
                return 0;
            }
            n_bodies = tmp;
        } catch (...) {
            std::cout << "Invalid argument for <number of bodies>" << std::endl;
        }
    }

    // first, initialize magick++
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
    for (size_t i = 0; i < n_bodies - 1; ++i) {
        bodies.m.push_back(10 + 5 * uniform());
        bodies.colors.push_back(colors[rand() % colors.size()]);
        bodies.r.emplace_back((0.25 + 0.5 * uniform()) * canvas_width,
                              (0.25 + 0.5 * uniform()) * canvas_height);

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
    switch (algo) {
        case 0:
            single_thread(bodies, n_threads, drawer);
            break;
        case 1:
            multi_thread_1(bodies, n_threads, drawer);
            break;
        case 2:
            multi_thread_2(bodies, n_threads, drawer);
            break;
        case 3:
            barnes_hut(bodies, n_threads, drawer);
        case 4:
            barnes_hut_multi(bodies, n_threads, drawer);
    }
#else
    switch (algo) {
        case 0:
            single_thread(bodies, n_threads);
            break;
        case 1:
            multi_thread_1(bodies, n_threads);
            break;
        case 2:
            multi_thread_2(bodies, n_threads);
            break;
        case 3:
            barnes_hut(bodies, n_threads);
        case 4:
            barnes_hut_multi(bodies, n_threads);
    }
#endif
    auto end = std::chrono::steady_clock::now();

    std::chrono::duration<double> duration = end - start;
    std::cout << "Simulation ended. Took " << duration.count() << "s"
              << std::endl;
    std::cout << "Average time per step: "
              << duration.count() / std::floor(t_end / dt) * 1000 << "ms"
              << std::endl;
}
