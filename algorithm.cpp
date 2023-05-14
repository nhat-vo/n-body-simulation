#include "Magick++/Color.h"
#include "Magick++/Drawable.h"
#include "Magick++/Image.h"
#include "Magick++/STL.h"
#include <Magick++.h>
#include <cmath>
#include <condition_variable>
#include <format>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

const double G = 6.67e-11; // gravitational constant

class Vect {
  public:
    double x, y;
    Vect(double x, double y) : x(x), y(y) {}
    double norm2() { return x * x + y * y; }
    double norm() { return sqrt(x * x + y * y); }
};
Vect operator+(const Vect &a, const Vect &b) { return {a.x + b.x, a.y + b.y}; }
Vect operator-(const Vect &a, const Vect &b) { return {a.x - b.x, a.y - b.y}; }
Vect operator*(const Vect &a, const Vect &b) { return {a.x * b.x, a.y * b.y}; }
Vect operator/(const Vect &a, const Vect &b) { return {a.x / b.x, a.y / b.y}; }

Vect operator+(const Vect &a, double b) { return {a.x + b, a.y + b}; }
Vect operator-(const Vect &a, double b) { return {a.x - b, a.y - b}; }
Vect operator*(const Vect &a, double b) { return {a.x * b, a.y * b}; }
Vect operator/(const Vect &a, double b) { return {a.x / b, a.y / b}; }

Vect operator+(double b, const Vect &a) { return {a.x + b, a.y + b}; }
Vect operator-(double b, const Vect &a) { return {a.x - b, a.y - b}; }
Vect operator*(double b, const Vect &a) { return {a.x * b, a.y * b}; }
Vect operator/(double b, const Vect &a) { return {a.x / b, a.y / b}; }

std::ostream &operator<<(std::ostream &os, const Vect &vect) {
    return os << "{" << vect.x << ", " << vect.y << "}";
}

class Body {
  public:
    const double m; // mass
    Vect r, v;      // position, velocity
    std::mutex lock;
    std::string color;

    Body(double mass, const Vect &initial_position,
         const Vect &initial_velocity, const std::string &color = "blue")
        : m(mass), r(initial_position), v(initial_velocity), color(color) {}
};

std::ostream &operator<<(std::ostream &os, const Body &body) {
    return os << "Body(" << body.m << ", " << body.r << ", " << body.v << ")";
}

void compute_forces(std::vector<Body *> &bodies, size_t start, size_t end,
                    double dt) {
    size_t n = bodies.size();
    for (size_t i = start; i < end; i++) {
        for (size_t j = i + 1; j < n; j++) {
            Vect dr = bodies[j]->r - bodies[i]->r;

            double dist_sq = dr.norm2();
            double dist = sqrt(dist_sq);
            double force_mag = G * bodies[i]->m * bodies[j]->m / dist_sq;
            Vect force = force_mag / dist * dr;

            bodies[i]->lock.lock();
            bodies[i]->v = bodies[i]->v + force / bodies[i]->m * dt;
            bodies[i]->lock.unlock();

            bodies[j]->lock.lock();
            bodies[j]->v = bodies[j]->v - force / bodies[j]->m * dt;
            bodies[j]->lock.unlock();
        }
    }
}

void update_positions(std::vector<Body *> &bodies, int start, int end,
                      double dt) {
    for (int i = start; i < end; i++)
        bodies[i]->r = bodies[i]->r + bodies[i]->v * dt;
}

void simulate(std::vector<Body *> &bodies, int num_threads, double dt,
              int num_steps) {
    int n = bodies.size();
    std::vector<std::thread> threads(num_threads);
    for (int step = 0; step < num_steps; step++) {
        // Compute forces
        int chunk_size = n / num_threads;
        int start = 0;
        int end = chunk_size;
        for (int i = 0; i < num_threads; i++) {
            if (i == num_threads - 1) {
                end = n;
            }
            threads[i] =
                std::thread(compute_forces, ref(bodies), start, end, dt);
            start = end;
            end += chunk_size;
        }
        for (auto &t : threads) {
            t.join();
        }

        // Update positions
        start = 0;
        end = chunk_size;
        for (int i = 0; i < num_threads; i++) {
            if (i == num_threads - 1) {
                end = n;
            }
            threads[i] =
                std::thread(update_positions, ref(bodies), start, end, dt);
            start = end;
            end += chunk_size;
        }
        for (auto &t : threads) {
            t.join();
        }
    }
}

std::mutex draw_lock;
std::condition_variable draw_signal;
const Vect canvas_size{800, 500};

void draw(std::vector<Magick::Image> &images, const std::vector<Body *> &bodies,
          double &draw_t, double &current_t, double end_t) {
    std::unique_lock<std::mutex> lock(draw_lock);
    while (draw_t < end_t) {
        while (draw_t < current_t)
            draw_signal.wait(lock);
        Magick::Image frame(Magick::Geometry(canvas_size.x, canvas_size.y),
                            Magick::Color("white"));
        for (const Body *body : bodies) {
            frame.fillColor(body->color);
            frame.draw(Magick::DrawableCircle(body->r.x, body->r.y,
                                              body->r.x - 10, body->r.y));
        }
        images.push_back(frame);
    }
}

int main(int argc, char **argv) {
    Magick::InitializeMagick(*argv);
    int num_bodies = 2;  // to be precised
    int num_threads = 2; // to be precised
    double dt = 1e-3;
    double t_end = 100;
    int start = 0;
    int end = 2;

    // double vy = sqrt(G * 1e14 / 100);
    const Vect offset = canvas_size / 2;
    double vy = 0;
    Body b1 = Body(1e14, offset + Vect{0, 0}, {0, vy}, "blue");
    Body b2 = Body(1e14, offset + Vect{100, 0}, {0, -vy}, "red");
    Body b3 =
        Body(1, offset + Vect{0, -50}, {sqrt(1.2 * G * 1e14 / 50), 0}, "green");
    std::cout << "body created" << std::endl;

    std::vector<Body *> bodies{&b1, &b3};

    std::vector<Magick::Image> images;
    double gif_dt = 1;
    double gif_t = 0;

    double t = 0;

    std::thread draw_thread(draw, std::ref(images), std::ref(bodies),
                            std::ref(gif_t), std::ref(t), t_end);

    for (; t <= t_end; t += dt) {
        compute_forces(bodies, start, end, dt);
        update_positions(bodies, start, end, dt);
        if (t >= gif_t) {
            std::lock_guard<std::mutex> _lk(draw_lock);
            gif_t += gif_dt;
            draw_signal.notify_one();
        }
    }
    draw_thread.join();
    Magick::writeImages(images.begin(), images.end(), "image.gif");
}
