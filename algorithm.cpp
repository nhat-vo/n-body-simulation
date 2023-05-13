#include <iostream>
#include <thread>
#include <vector>
#include <cmath>
#include <mutex>

using namespace std;

const double G = 6.67e-11; // gravitational constant

struct Body {
    double mass;
    double x, y;
    double vx, vy;
};

void compute_forces(vector<Body>& bodies, int start, int end, vector<mutex>& mutexes) {
    int n = bodies.size();
    for (int i = start; i < end; i++) {
        for (int j = i + 1; j < n; j++) {
            double dx = bodies[j].x - bodies[i].x;
            double dy = bodies[j].y - bodies[i].y;
            double dist_sq = dx * dx + dy * dy;
            double dist = sqrt(dist_sq);
            double force = G * bodies[i].mass * bodies[j].mass / dist_sq;
            double fx = force * dx / dist;
            double fy = force * dy / dist;
            mutexes[i].lock();
            bodies[i].vx += fx / bodies[i].mass;
            bodies[i].vy += fy / bodies[i].mass;
            mutexes[i].unlock();
            mutexes[j].lock();
            bodies[j].vx -= fx / bodies[j].mass;
            bodies[j].vy -= fy / bodies[j].mass;
            mutexes[j].unlock();
        }
    }
}

void update_positions(vector<Body>& bodies, int start, int end, double dt) {
    for (int i = start; i < end; i++) {
        bodies[i].x += bodies[i].vx * dt;
        bodies[i].y += bodies[i].vy * dt;
    }
}

void simulate(vector<Body>& bodies, int num_threads, double dt, int num_steps) {
    int n = bodies.size();
    vector<thread> threads(num_threads);
    vector<mutex> mutexes(n);
    for (int step = 0; step < num_steps; step++) {
        // Compute forces
        int chunk_size = n / num_threads;
        int start = 0;
        int end = chunk_size;
        for (int i = 0; i < num_threads; i++) {
            if (i == num_threads - 1) {
                end = n;
            }
            threads[i] = thread(compute_forces, ref(bodies), start, end, ref(mutexes));
            start = end;
            end += chunk_size;
        }
        for (auto& t : threads) {
            t.join();
        }

        // Update positions
        start = 0;
        end = chunk_size;
        for (int i = 0; i < num_threads; i++) {
            if (i == num_threads - 1) {
                end = n;
            }
            threads[i] = thread(update_positions, ref(bodies), start, end, dt);
            start = end;
            end += chunk_size;
        }
        for (auto& t : threads) {
            t.join();
        }
    }
}

int main() {
    int num_bodies = 2; // to be precised
    int num_threads = 2; // to be precised
    double dt = 1.0;
    int num_steps = 10; 
    int start = 0;
    int end = 10;

    vector<Body> bodies(num_bodies);
    bodies[0].mass = 1;
    bodies[0].x = 2;
    bodies[1].mass = 2;
    bodies[1].x = 2;

    for (int step = 0; step < num_steps; ++step) {
        compute_forces(bodies, start, end, mutexes);
        update_positions(bodies, start, end, dt);
    }
}




