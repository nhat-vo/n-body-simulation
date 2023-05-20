#include <Magick++.h>
#include <atomic>
#include <barrier>
#include <cmath>
#include <condition_variable>
#include <format>
#include <functional>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

const double G = 6.67e-11; // gravitational constant
const double dt = 1e-3;
const double t_end = 100;
const double draw_dt = 1;
const size_t canvas_width = 800, canvas_height = 500;
const size_t n_threads = 2;

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
            for (size_t i = 0; i < drawer.draw_r.size(); ++i) {
                const Vect &r = drawer.draw_r[i];
                drawer.images.back().fillColor(drawer.colors[i]);
                drawer.images.back().draw(
                    Magick::DrawableCircle(r.x, r.y, r.x - 5, r.y));
            }
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
    std::vector<Vect> v;
    std::vector<std::string> colors;
};

void compute_forces(Scenario &bodies, size_t start, size_t end, double dt) {
    size_t n = bodies.r.size();
    for (size_t i = start; i < end; i++) {
        for (size_t j = i + 1; j < n; j++) {
            Vector dr = bodies.r[j] - bodies.r[i];

            double dist_sq = dr.norm2();
            double dist = sqrt(dist_sq);
            double force_mag = G * bodies.m[i] * bodies.m[j] / dist_sq;
            Vector force = dr * (force_mag / dist);

            bodies.v[i] = bodies.v[i] + force * (dt / bodies.m[i]);
            bodies.v[j] = bodies.v[j] - force * (dt / bodies.m[j]);
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

// void multi_thread_aux(ThreadArgs &args, size_t start, size_t end) {
//     size_t n = args.curr_r->size();
//     for (double t = 0; t <= t_end; t += dt) {
//         for (size_t i = start; i < end; i++) {
//             for (size_t j = i + 1; j < n; j++) {
//                 Vector dr = (*args.curr_r)[j] - (*args.curr_r)[i];

//                 double dist_sq = dr.norm2();
//                 double dist = sqrt(dist_sq);
//                 double force_mag = G * args.m[i] * args.m[j] / dist_sq;
//                 Vector force = dr * (force_mag / dist);

//                 args.v[i] += force * (dt * dt / args.m[i]);
//                 args.v[j] -= force * (dt * dt / args.m[j]);
//             }
//         }

//         args.barrier.arrive_and_wait();
//         std::cout << t << std::endl;

//         for (size_t i = start; i < end; i++) {
//             (*args.next_r)[i] = (*args.curr_r)[i] + args.v[i];
//         }

//         args.barrier.arrive_and_wait();
//         std::cout << t << std::endl;
//     }
// }

// void multi_thread(Scenario &bodies, double &draw_t, std::vector<Vect>
// &draw_pos,
//                   bool &draw_pending) {
//     size_t n = bodies.r.size();
//     std::vector<std::thread> threads(n_threads - 1);
//     std::vector<Vect> r1 = bodies.r, r2(n);
//     std::vector<Vect> *curr_r = &r1, *next_r = &r2;
//     std::vector<AVect> v(n);
//     for (size_t i = 0; i < n; ++i) {
//         v[i].x = bodies.v[i].x;
//         v[i].y = bodies.v[i].y;
//     }

//     size_t chunk_size = n / n_threads;
//     std::barrier barrier(n_threads);
//     ThreadArgs args{.curr_r = curr_r,
//                     .next_r = next_r,
//                     .v = v,
//                     .m = bodies.m,
//                     .common_t = 0,
//                     .pending_count = n_threads,
//                     .draw_pending = draw_pending,
//                     .draw_lock = draw_lock,
//                     .draw_signal = draw_signal,
//                     .draw_pos = draw_pos,
//                     .draw_t = draw_t,
//                     .barrier = barrier};
//     for (size_t i = 0; i < n_threads - 1; ++i) {
//         threads[i] = std::thread(multi_thread_aux, std::ref(args),
//                                  i * chunk_size, (i + 1) * chunk_size);
//     }
//     multi_thread_aux(args, (n_threads - 1) * chunk_size, n);

//     for (size_t i = 0; i < n_threads - 1; ++i) {
//         threads[i].join();
//     }
// }

int main(int argc, char **argv) {
    Magick::InitializeMagick(*argv);

    // setup scenario
    const Vect offset = {canvas_width / 2, canvas_height / 2};
    std::cout << "body created" << std::endl;
    Scenario bodies{
        {1e14, 1},
        {offset + Vect(0, 0), offset + Vect(0, -50)},
        {{0, 0}, {sqrt(1.2 * G * 1e14 / 50), 0}},
        {"red", "green"},
    };

    // setup drawing thread
    Drawer drawer(bodies.colors);

    single_thread(bodies, drawer);

    // drawer.terminate();
}
