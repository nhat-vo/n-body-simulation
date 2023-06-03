#include "algorithms.hpp"
static void compute_forces(Scenario &bodies, size_t start, size_t end,
                           double dt) {
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

static void update_positions(Scenario &bodies, int start, int end, double dt) {
    for (int i = start; i < end; i++)
        bodies.r[i] = bodies.r[i] + bodies.v[i] * dt;
}

#ifdef VISUALIZE
void single_thread(Scenario &bodies, size_t n_threads, Drawer &drawer) {
#else
void single_thread(Scenario &bodies, size_t n_threads) {
#endif

    std::cout << "Simulation with single_thread\n";

    for (double t = 0; t < t_end; t += dt) {
#ifdef VISUALIZE
        drawer.trigger_draw(t, &bodies.r);
#endif

        compute_forces(bodies, 0, bodies.r.size(), dt);
        update_positions(bodies, 0, bodies.r.size(), dt);
    }
}
