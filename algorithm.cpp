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
