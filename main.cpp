#include <chrono>
#include <cmath>

#include "headers/algorithms.hpp"

int main(int argc, char** argv) {
    std::string algos[5] = {"single-thread", "multi-thread-1", "multi-thread-2",
                            "barnes-hut", "barnes-hut multi-threaded"};

    // nothing fancy here, just parsing command line arguments
    if (argc == 1) {
        std::cout << "Usage: " << argv[0] << " <algorithm id (default = 0)>"
                  << " <number of threads (default = 1)>"
                  << " <number of bodies (default = 100)>" << std::endl;

        std::cout << "Algorithms:" << std::endl;
        for (int i = 0; i < 5; i++) {
            std::cout << i << ": " << algos[i] << std::endl;
        }

        return 0;
    }

    if (argc > 1) {
        // Iterate over the arguments
        for (int i = 1; i < argc; ++i) {
            // Compare the argument with the desired value
            if (std::string(argv[i]) == "--benchmark") {
                std::cout << "Found desired_value in argv!" << std::endl;
                // run benchmark
                return 0;
            }
        }
    }

    size_t n_threads = 1, n_bodies = 100, algo = 4;
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

#ifdef DEBUG
    std::cout << "Simulating " << n_bodies << " bodies with " << n_threads
              << " threads for " << floor(t_end / dt) << " steps." << std::endl;
#else
    std::cout << algos[algo] << "; Bodies:" << n_bodies
              << "; Threads:" << n_threads;
#endif

#ifdef WRITE
    Writer writer_single("single.txt");
    Writer writer_barnes("barnes.txt");
#endif

    // setup scenario
    const Vect offset = {canvas_width / 2, canvas_height / 2};

    Scenario bodies{{1e14}, {offset + Vect(0, 0)}, {{0, 0}}, {"red"}};
    std::vector<std::string> colors{"blue", "green",  "gold",   "grey",
                                    "pink", "orange", "purple", "brown"};
    initialize_bodies(bodies, n_bodies, colors);

#ifdef VISUALIZE
    // Drawing setup
    Magick::InitializeMagick(*argv);
    Drawer drawer(bodies.colors);
#endif

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
            break;
        case 4:
            barnes_hut_multi(bodies, n_threads, drawer);
            break;
    }
#elif WRITE
    Scenario bodies2 = bodies;
    single_thread(bodies, n_threads, writer_single);
    barnes_hut(bodies2, n_threads, writer_barnes);
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
            break;
        case 4:
            barnes_hut_multi(bodies, n_threads);
            break;
    }
#endif
    auto end = std::chrono::steady_clock::now();

    std::chrono::duration<double> duration = end - start;

#ifdef DEBUG
    std::cout << "Time: " << duration.count() << "s" << std::endl;
    std::cout << "Average time per step: "
              << duration.count() / std::floor(t_end / dt) * 1000 << "ms"
              << std::endl;
#else
    std::cout << "; Time: " << duration.count() << "s" << std::endl;
#endif
}

void initialize_bodies(Scenario& bodies, size_t n_bodies,
                       std::vector<std::string>& colors) {
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
}
