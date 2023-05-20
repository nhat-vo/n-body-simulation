#include <Magick++.h>
#include <atomic>
#include <barrier>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <ctime>
#include <format>
#include <functional>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

const double G = 6.67e-11;
const double PI = 3.1415926535;
const size_t N_BODIES = 100;

inline double uniform() { return (double)rand() / RAND_MAX; }
inline size_t discrete_uniform(size_t n) { return rand() % n; }

namespace config {
const double dt = 1e-3;
const double t_end = 100;
const double draw_dt = 1;
const size_t canvas_width = 400, canvas_height = 400;
const size_t n_threads = 2;
} // namespace config

using namespace config;

template <typename T> class Vector {
  public:
    double x, y;
    Vector(double x, double y) : x(x), y(y) {}
    Vector() : x(0), y(0) {}
    double norm2() { return x * x + y * y; }
    double norm() { return sqrt(x * x + y * y); }

    template <typename U> void operator+=(const Vector<U> &b) {
        x += b.x;
        y += b.y;
    }
    template <typename U> void operator-=(const Vector<U> &b) {
        x -= b.x;
        y -= b.y;
    }
    template <typename U> Vector operator+(const Vector<U> &b) const {
        return {x + b.x, y + b.y};
    }
    template <typename U> Vector operator-(const Vector<U> &b) const {
        return {x - b.x, y - b.y};
    }
    template <typename U> Vector operator*(const Vector<U> &b) const {
        return {x * b.x, y * b.y};
    }
    template <typename U> Vector operator/(const Vector<U> &b) const {
        return {x / b.x, y / b.y};
    }

    Vector operator+(double b) { return {x + b, y + b}; }
    Vector operator-(double b) { return {x - b, y - b}; }
    Vector operator*(double b) { return {x * b, y * b}; }
    Vector operator/(double b) { return {x / b, y / b}; }
};

typedef Vector<double> Vect;
typedef Vector<std::atomic<double>> AVect;

class Drawer {
  public:
    std::vector<Magick::Image> images;
    double draw_t = 0;
    bool draw_pending = false;
    std::vector<Vect> draw_r;
    const std::vector<std::string> &colors;
    std::thread draw_thread;
    std::mutex draw_lock;
    std::condition_variable draw_signal;

    Drawer(const std::vector<std::string> &colors)
        : colors(colors), draw_thread(&Drawer::draw, this) {}

    static void draw(Drawer *drawer_pt) {
        Drawer &drawer = *drawer_pt;

        std::unique_lock<std::mutex> lock(drawer.draw_lock);
        const Magick::Geometry image_size(canvas_width, canvas_height);
        const Magick::Color background_color("white");

        while (drawer.draw_t < t_end) {
            while (!drawer.draw_pending)
                drawer.draw_signal.wait(lock);

            drawer.images.emplace_back(image_size, background_color);
            auto &frame = drawer.images.back();

            for (size_t i = 0; i < drawer.draw_r.size(); ++i) {
                const Vect &r = drawer.draw_r[i];
                double x = r.x, y = (double)canvas_height - r.y;
                frame.fillColor(drawer.colors[i]);
                frame.draw(Magick::DrawableCircle(x, y, x - 5, y));
            }

            frame.strokeColor("black");
            frame.draw(Magick::DrawableText(
                10, canvas_height - 10,
                "t = " + std::to_string((int)drawer.draw_t)));

            drawer.draw_pending = false;
            drawer.draw_signal.notify_all();
        }
    }
    void trigger_draw(double t, std::vector<Vect> *curr_r) {
        if (t >= draw_t) {
            std::unique_lock<std::mutex> lock(draw_lock);

            // wait until the draw thread has finished drawing
            while (draw_pending)
                draw_signal.wait(lock);
            draw_t += draw_dt;
            draw_r = *curr_r;

            // trigger the drawing thread
            draw_pending = true;
            draw_signal.notify_one();
        }
    }
    ~Drawer() {
        draw_thread.join();
        std::cout << "rendering done, exporting images" << std::endl;
        Magick::writeImages(images.begin(), images.end(), "image.gif");
        std::cout << "image written to "
                  << "image.gif" << std::endl;
    }
};

struct Scenario {
    std::vector<double> m;
    std::vector<Vect> r;
    std::vector<AVect> v;
    std::vector<std::string> colors;
};

void compute_forces(Scenario &bodies, size_t start, size_t end, double dt) {
    size_t n = bodies.r.size();
    for (size_t i = start; i < end; i++) {
        for (size_t j = i + 1; j < n; j++) {
            Vector dr = bodies.r[j] - bodies.r[i];

            double dist_sq = std::max(dr.norm2(), 1e-6);
            double dist = sqrt(dist_sq);
            double force_mag = G * bodies.m[i] * bodies.m[j] / dist_sq;
            Vector force = dr * (force_mag / dist);

            bodies.v[i] += force * (dt / bodies.m[i]);
            bodies.v[j] -= force * (dt / bodies.m[j]);
        }
    }
}

void update_positions(Scenario &bodies, int start, int end, double dt) {
    for (int i = start; i < end; i++)
        bodies.r[i] = bodies.r[i] + bodies.v[i] * dt;
}

void single_thread(Scenario &bodies, Drawer &drawer) {
    for (double t = 0; t <= t_end; t += dt) {
        drawer.trigger_draw(t, &bodies.r);

        compute_forces(bodies, 0, bodies.r.size(), dt);
        update_positions(bodies, 0, bodies.r.size(), dt);
    }
}

void multi_thread_naive(Scenario &bodies, Drawer &drawer) {
    size_t n = bodies.r.size();
    size_t chunk_size = n / n_threads;
    if (chunk_size == 0) {
        single_thread(bodies, drawer);
        return;
    }
    std::vector<std::thread> threads(n_threads - 1);
    std::vector<Vect> r1 = bodies.r, r2(n);
    std::vector<AVect> v(n);
    for (size_t i = 0; i < n; ++i) {
        v[i].x = bodies.v[i].x;
        v[i].y = bodies.v[i].y;
    }

    for (double t = 0; t <= t_end; t += dt) {
        drawer.trigger_draw(t, &bodies.r);
        for (size_t i = 0; i < n_threads - 1; ++i) {
            threads[i] = std::thread(compute_forces, std::ref(bodies),
                                     i * chunk_size, (i + 1) * chunk_size, dt);
        }
        compute_forces(bodies, (n_threads - 1) * chunk_size, n, dt);
        for (size_t i = 0; i < n_threads - 1; ++i)
            threads[i].join();

        for (size_t i = 0; i < n_threads - 1; ++i) {
            threads[i] = std::thread(update_positions, std::ref(bodies),
                                     i * chunk_size, (i + 1) * chunk_size, dt);
        }
        update_positions(bodies, (n_threads - 1) * chunk_size, n, dt);
        for (size_t i = 0; i < n_threads - 1; ++i)
            threads[i].join();
    }
}

int main(int argc, char **argv) {
    Magick::InitializeMagick(*argv);

    // setup scenario
    const Vect offset = {canvas_width / 2, canvas_height / 2};
    std::cout << "body created" << std::endl;

    double v = sqrt(1 * G * 1e14 / 50);
    // Scenario bodies{
    //     {1e14, 1},
    //     {offset + Vect(0, 0), offset + Vect(0, -50)},
    //     {{0, 0}, {v, 0}},
    //     {"red", "green"},
    // };
    Scenario bodies{{1e14}, {offset + Vect(0, 0)}, {{0, 0}}, {"red"}};
    std::vector<std::string> colors{"blue", "green",  "gold",   "grey",
                                    "pink", "orange", "purple", "brown"};
    for (size_t i = 0; i < N_BODIES; ++i) {
        bodies.m.push_back(10 + 5 * uniform());
        bodies.colors.push_back(colors[rand() % colors.size()]);
        bodies.r.push_back(Vect((0.5 + 0.2 * uniform()) * canvas_width,
                                (0.5 + 0.2 * uniform()) * canvas_height));
        Vect dir = bodies.r.back() - bodies.r.front();
        double dist = dir.norm();
        double v = sqrt((0.5 + uniform()) * G * bodies.m.front() / dist);
        dir = dir / dist;
        bodies.v.push_back({v * (-dir.y), v * dir.x});
    }

    // setup drawing thread
    Drawer drawer(bodies.colors);

    std::cout << "Staring simulation..." << std::endl;
    auto start = std::chrono::steady_clock::now();
    multi_thread_naive(bodies, drawer);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Simulation ended. Took " << duration.count() << "s"
              << std::endl;
}
