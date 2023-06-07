#include "barnes-hut.hpp"
#include <vector>

#ifdef VISUALIZE
void barnes_hut_multi(Scenario &bodies, size_t n_threads, Drawer &drawer) {
#elif WRITE
void barnes_hut_multi(Scenario &bodies, size_t n_threads, Writer &writer) {
#else
void barnes_hut_multi(Scenario &bodies, size_t n_threads) {
#endif

    std::cout << "Simulation with multi-threaded Barnes-Hut\n";

    for (double t = 0; t < t_end; t += dt) {
#ifdef VISUALIZE
        drawer.trigger_draw(t, &bodies.r);
#elif WRITE
        writer.write(bodies.r);
#endif

        barnes_hut_update_step_multi(bodies, n_threads);
    }
}

void barnes_hut_update_step_aux(int start, int end, Scenario &bodies, QuadNode *root) {
    for (int i = start; i < end; i++) {
        const double m = bodies.m[i];
        const Vect &r = bodies.r[i];

        auto update_v = [&](const Vect &other_r, const double other_m) {
            Vect dr = other_r - r;
            double dist_sq = std::max(dr.norm2(), 1e-6);
            double dist = sqrt(dist_sq);
            double force_mag = G * other_m * m / dist_sq;
            Vect force = dr * (force_mag / dist);
            bodies.v[i] += force * (dt / m);
        };

        std::stack<QuadNode *> stack;
        stack.push(root);
        while (!stack.empty()) {
            QuadNode *curr = stack.top();
            stack.pop();

            if (!curr->body_id.empty()) {
                for (int curr_body : curr->body_id) {
                    if (curr_body != i) {
                        update_v(bodies.r[curr_body], bodies.m[curr_body]);
                    }
                }
            } else if (curr->isFarEnough(bodies.r[i])) {
                update_v(curr->center_of_mass, curr->m);
            } else {
                for (int j = 0; j < 4; j++) {
                    if (curr->children[j]) stack.push(curr->children[j]);
                }
            }
        }
    }

}

void barnes_hut_update_step_multi(Scenario &bodies, int num_threads) {
    QuadNode *root = QuadNode::constructBarnesHutTree(&bodies);
    std::thread threads[num_threads-1];
    int num_bodies = bodies.r.size();
    int num_bodies_per_thread = num_bodies / num_threads;

    for (int i=0; i < num_threads-1; i++)
        threads[i] = std::thread(barnes_hut_update_step_aux, i*num_bodies_per_thread, (i+1)*num_bodies_per_thread, std::ref(bodies), root);
    barnes_hut_update_step_aux((num_threads-1)*num_bodies_per_thread, num_bodies, bodies, root);

    for (auto &thread : threads)
        thread.join();

    /* Update positions */
    for (size_t i = 0; i < bodies.r.size(); i++) {
        bodies.r[i] += bodies.v[i] * dt;
    }
    delete root;
}