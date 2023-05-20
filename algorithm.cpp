// #define VISUALIZE
#ifdef VISUALIZE
#include <Magick++.h>
#endif

#include <chrono>
#include <cmath>
#include <iostream>
#include <list>
#include <thread>
#include <vector>

const double G = 6.67e-11;
const double PI = 3.1415926535;
const size_t N_BODIES = 1000;

inline double uniform() { return (double)rand() / RAND_MAX; }
inline size_t discrete_uniform(size_t n) { return rand() % n; }

namespace config {
const double dt = 1e-2;
const double t_end = 100;
const double draw_dt = 1;
const size_t canvas_width = 400, canvas_height = 400;
const size_t n_threads = 1;
} // namespace config

using namespace config;

class Vector {
  public:
    double x, y;
    Vector(double x, double y) : x(x), y(y) {}
    Vector() : x(0), y(0) {}
    double norm2() { return x * x + y * y; }
    double norm() { return sqrt(x * x + y * y); }

    void operator+=(const Vector &b) {
        x += b.x;
        y += b.y;
    }
    void operator-=(const Vector &b) {
        x -= b.x;
        y -= b.y;
    }
    Vector operator+(const Vector &b) const { return {x + b.x, y + b.y}; }
    Vector operator-(const Vector &b) const { return {x - b.x, y - b.y}; }
    Vector operator*(const Vector &b) const { return {x * b.x, y * b.y}; }
    Vector operator/(const Vector &b) const { return {x / b.x, y / b.y}; }

    Vector operator+(double b) { return {x + b, y + b}; }
    Vector operator-(double b) { return {x - b, y - b}; }
    Vector operator*(double b) { return {x * b, y * b}; }
    Vector operator/(double b) { return {x / b, y / b}; }
};

typedef Vector Vect;

#ifdef VISUALIZE
class Drawer {
  private:
    const Magick::Geometry image_size;
    const Magick::Color background_color;
    const std::vector<std::string> &colors;

    std::list<Magick::Image> images;
    double draw_t = 0;
    std::vector<Vect> draw_r;

    std::thread draw_thread;

    static void draw(Drawer *drawer_pt) {
        // auto start = std::chrono::steady_clock::now();
        Drawer &drawer = *drawer_pt;

        drawer.images.emplace_back(drawer.image_size, drawer.background_color);
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
            "t = " + std::to_string((int)(drawer.draw_t - draw_dt))));
        // auto end = std::chrono::steady_clock::now();
        // std::chrono::duration<double> duration = end - start;
        // std::cout << "draw time: " << duration.count() << std::endl;
    }

  public:
    Drawer(const std::vector<std::string> &colors)
        : image_size(canvas_width, canvas_height), background_color("white"),
          colors(colors) {}

    ~Drawer() {
        if (draw_thread.joinable())
            draw_thread.join();
        std::cout << "rendering done, exporting images" << std::endl;
        Magick::writeImages(images.begin(), images.end(), "image.gif");
        std::cout << "image written to "
                  << "image.gif" << std::endl;
    }

    void trigger_draw(double t, std::vector<Vect> *curr_r) {
        if (t >= draw_t) {
            // wait until the draw thread has finished drawing
            if (draw_thread.joinable())
                draw_thread.join();

            draw_t += draw_dt;
            draw_r = *curr_r;

            draw_thread = std::thread(Drawer::draw, this);
        }
    }
};
#endif

struct Scenario {
    std::vector<double> m;
    std::vector<Vect> r;
    std::vector<Vect> v;
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

#ifdef VISUALIZE
void single_thread(Scenario &bodies, Drawer &drawer) {
#else
void single_thread(Scenario &bodies) {
#endif
    for (double t = 0; t <= t_end; t += dt) {
#ifdef VISUALIZE
        drawer.trigger_draw(t, &bodies.r);
#endif

        compute_forces(bodies, 0, bodies.r.size(), dt);
        update_positions(bodies, 0, bodies.r.size(), dt);
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
void multi_thread_1(Scenario &bodies, Drawer &drawer) {
#else
void multi_thread_1(Scenario &bodies) {
#endif
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
#ifdef VISUALIZE
    Magick::InitializeMagick(*argv);
#endif

    // setup scenario
    const Vect offset = {canvas_width / 2, canvas_height / 2};
    std::cout << "body created" << std::endl;

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
    for (size_t i = 0; i < N_BODIES; ++i) {
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

    std::cout << "Staring simulation..." << std::endl;

    auto start = std::chrono::steady_clock::now();
#ifdef VISUALIZE
    multi_thread_1(bodies, drawer);
#else
    multi_thread_1(bodies);
#endif
    auto end = std::chrono::steady_clock::now();

    std::chrono::duration<double> duration = end - start;
    std::cout << "Simulation ended. Took " << duration.count() << "s"
              << std::endl;
}
