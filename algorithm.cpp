#include <cmath>
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

    Body(double mass, const Vect &initial_position,
         const Vect &initial_velocity)
        : m(mass), r(initial_position), v(initial_velocity) {}
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

            // bodies[i].lock.lock();
            bodies[i]->v = bodies[i]->v + force / bodies[i]->m * dt;
            // bodies[i].lock.unlock();

            // bodies[j].lock.lock();
            bodies[j]->v = bodies[j]->v - force / bodies[j]->m * dt;
            // bodies[j].lock.unlock();
        }
    }
}

void update_positions(std::vector<Body *> &bodies, int start, int end,
                      double dt) {
    for (int i = start; i < end; i++)
        bodies[i]->r = bodies[i]->r + bodies[i]->v * dt;
}

void simulate(std::vector<Body> &bodies, int num_threads, double dt,
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

int main() {
    int num_bodies = 2;  // to be precised
    int num_threads = 2; // to be precised
    double dt = 1e-3;
    double t_end = 10;
    int start = 0;
    int end = 2;

    Body b1 = Body(1e13, {-10, 0}, {0, 0});
    Body b2 = Body(1e13, {10, 0}, {0, 0});
    std::cout << "body created" << std::endl;

    std::vector<Body *> bodies{&b1, &b2};

    for (double t = 0; t < t_end; t += dt) {
        compute_forces(bodies, start, end, dt);
        update_positions(bodies, start, end, dt);
        std::cout << t << ": " << *bodies[0] << " " << *bodies[1] << std::endl;
        if (bodies[0]->r.x >= 0)
            break;
    }
}
