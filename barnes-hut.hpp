#include <chrono>
#include <cmath>
#include <iostream>
#include <stack>
#include <thread>
#include <tuple>
#include <vector>

#include "common.hpp"
#include "vect.hpp"
#include "visualizer.hpp"

using namespace config;

class QuadNode {
    enum quad { nw, ne, sw, se };  // indeed this is not used
    bool is_empty = true;

public:
    QuadNode *children[4]{nullptr, nullptr, nullptr, nullptr};
    double m = 0;
    Vect center_of_mass;
    const Vect center;  // center and dimension should be const, right?
    const Vect dimension;
    Scenario *const scenario;
    std::vector<int> body_id;

    // This is the main entry point of the Barnes-Hut tree. This constructs a
    // Barnes-Hut tree from `bodies`.
    // NOTE:: Remember to delete the result.
    static QuadNode *constructBarnesHutTree(Scenario *bodies) {
        QuadNode *root =
            new QuadNode(bodies, Vect(canvas_width / 2., canvas_height / 2.),
                         Vect(universe_width, universe_height));

        for (int i = 0; i < bodies->r.size(); i++) root->addBody(i);

        return root;
    }

    QuadNode(Scenario *const bodies, const Vect &center, const Vect &dimension)
        : center_of_mass(center),
          center(center),
          dimension(dimension),
          scenario(bodies) {}

    ~QuadNode() {
        for (int i = 0; i < 4; i++) delete children[i];
    }

    bool isFarEnough(Vect point) {
        Vect dr = center_of_mass - point;
        double dist_sq = std::max(dr.norm2(), 1e-6);
        double dist = sqrt(dist_sq);
        return dimension.x / dist < theta;
    }

private:
    quad getQuad(const Vect &r) const {
        return r.x < center.x ? (r.y < center.y ? sw : nw)
                              : (r.y < center.y ? se : ne);
    }

    void updateCenterOfMass(size_t id) {
        double new_m = m + scenario->m[id];
        center_of_mass =
            (center_of_mass * m + scenario->r[id] * scenario->m[id]) / new_m;
        m = new_m;
    }

    Vect getQuadCenter(quad quad) const {
        switch (quad) {
            case nw:
                return center + Vect(-dimension.x / 4, dimension.y / 4);
            case ne:
                return center + Vect(dimension.x / 4, dimension.y / 4);
            case sw:
                return center + Vect(-dimension.x / 4, -dimension.y / 4);
            case se:
                return center + Vect(dimension.x / 4, -dimension.y / 4);
            default:
                return Vect(0, 0);
        };
    }

    // Recursively add body `index` to this node
    void addBody(int index) {
        if (!isInside(scenario->r[index])) {
            std::cout << "Body " << index << " is not inside the node\n";
            std::cout << scenario->r[index] << " " << center << " " << dimension
                      << std::endl;
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
            children[quad_old] =
                new QuadNode(scenario, quad_center, dimension / 2);
            children[quad_old]->addBody(body_id[0]);
            body_id.clear();
        }

        // Get quad where new body belongs
        quad quad_new = getQuad(scenario->r[index]);

        // If the new node is empty, initialize the node
        if (children[quad_new] == nullptr) {
            Vect quad_center = getQuadCenter(quad_new);
            children[quad_new] =
                new QuadNode(scenario, quad_center, dimension / 2);
        }
        // Add the body to the node and update center of mass
        children[quad_new]->addBody(index);
        updateCenterOfMass(index);
    }

    bool isInside(const Vect &point) const {
        return point.x <= center.x + dimension.x / 2 &&
               point.x >= center.x - dimension.x / 2 &&
               point.y <= center.y + dimension.y / 2 &&
               point.y >= center.y - dimension.y / 2;
    }
};

void barnes_hut_update_step(Scenario &bodies) {
    QuadNode *root = QuadNode::constructBarnesHutTree(&bodies);

    /* Calculate the force exerted */
    for (int i = 0; i < bodies.r.size(); i++) {
        const double m = bodies.m[i];
        const Vect &r = bodies.r[i];

        auto update_v = [&bodies, m, &r, i](const Vect &other_r,
                                            const double other_m) {
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

    /* Update positions */
    for (size_t i = 0; i < bodies.r.size(); i++) {
        bodies.r[i] += bodies.v[i] * dt;
    }
    delete root;
}
