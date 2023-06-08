#include "../headers/barnes-hut.hpp"

#include "../headers/algorithms.hpp"
#ifdef VISUALIZE
void barnes_hut(Scenario &bodies, size_t n_threads, Drawer &drawer) {
#elif defined(WRITE)
void barnes_hut(Scenario &bodies, size_t n_threads, Writer &writer) {
#else
void barnes_hut(Scenario &bodies, size_t n_threads) {
#endif

#ifdef DEBUG
    std::cout << "Simulation with single-threaded Barnes-Hut\n";
#endif

    for (double t = 0; t < t_end; t += dt) {
#ifdef VISUALIZE
        drawer.trigger_draw(t, &bodies.r);
#elif WRITE
        writer.write(bodies.r);
#endif

        barnes_hut_update_step(bodies);
    }
}

void barnes_hut_update_step(Scenario &bodies) {
    QuadNode *root = QuadNode::constructBarnesHutTree(&bodies);

    /* Calculate the force exerted */
    for (size_t i = 0; i < bodies.r.size(); i++) {
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
                for (size_t curr_body : curr->body_id) {
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

    /* Update positions */
    for (size_t i = 0; i < bodies.r.size(); i++) {
        bodies.r[i] += bodies.v[i] * dt;
    }
    delete root;
}

void QuadNode::addBody(int index) {
    if (!isInside(scenario->r[index])) {
#ifdef DEBUG
        std::cout << "Body " << index << " is not inside the node\n";
        std::cout
            << "This is either because the body exitted the universe, or I "
               "did something wrong, so use the following to debug\n";
        std::cout << "Body position: " << scenario->r[index] << std::endl;
        std::cout << "Node center: " << center << std::endl;
        std::cout << "Node dimension: " << dimension << std::endl;
#endif
        return;
    }

    // If node is empty, add body to node
    if (is_empty || dimension.x < 1e-3) {
        is_empty = false;
        body_id.push_back(index);
        updateCenterOfMass(index);
        return;
    }

    // If there is only one body in the node, create a new define the sub
    // node of the old because the node was not split yet. The case where
    // there are multiple bodies is handled above, so there should always be
    // only one body here.
    if (!body_id.empty()) {
        quad quad_old = getQuad(scenario->r[body_id[0]]);
        Vect quad_center = getQuadCenter(quad_old);
        children[quad_old] = new QuadNode(scenario, quad_center, dimension / 2);
        children[quad_old]->addBody(body_id[0]);
        body_id.clear();
    }

    // Get quad where new body belongs
    quad quad_new = getQuad(scenario->r[index]);

    // If the new node is empty, initialize the node
    if (children[quad_new] == nullptr) {
        Vect quad_center = getQuadCenter(quad_new);
        children[quad_new] = new QuadNode(scenario, quad_center, dimension / 2);
    }
    // Add the body to the node and update center of mass
    children[quad_new]->addBody(index);
    updateCenterOfMass(index);
}
