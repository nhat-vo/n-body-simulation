#include <chrono>
#include <cmath>
#include <iostream>
#include <stack>
#include <thread>
#include <tuple>
#include <vector>
#include <omp.h>

#include "common.hpp"
#include "vect.hpp"
#include "visualizer.hpp"

using namespace config;

void barnes_hut_update_step(Scenario &bodies);
void barnes_hut_update_step_multi(Scenario &bodies, int num_threads);

class QuadNode {
    enum quad { nw, ne, sw, se };  // indeed this is not used
    bool is_empty = true;
    const Vect center;  // center and dimension should be const, right?
    const Vect dimension;
    Scenario *const scenario;

public:
    QuadNode *children[4]{nullptr, nullptr, nullptr, nullptr};
    double m = 0;
    Vect center_of_mass;
    std::vector<int> body_id;

    // This is the main entry point of the Barnes-Hut tree. This constructs a
    // Barnes-Hut tree from `bodies`.
    // NOTE:: Remember to delete the result.
    static QuadNode *constructBarnesHutTree(Scenario *bodies) {
        QuadNode *root =
            new QuadNode(bodies, Vect(canvas_width / 2., canvas_height / 2.),
                         Vect(universe_width, universe_height));

        #pragma omp parallel for
        for (size_t i = 0; i < bodies->r.size(); i++) {
            #pragma omp critical
            root->addBody(i);
        }

        return root;
    }

    QuadNode(Scenario *const bodies, const Vect &center, const Vect &dimension)
        : center(center),
          dimension(dimension),
          scenario(bodies),
          center_of_mass(center) {}

    ~QuadNode() {
        for (int i = 0; i < 4; i++) delete children[i];
    }

    bool isFarEnough(const Vect &point) const {
        Vect dr = center_of_mass - point;
        double dist_sq = std::max(dr.norm2(), 1e-6);
        return dimension.x * dimension.x / dist_sq < theta * theta;
    }

private:
    // Recursively add body `index` to this node
    void addBody(int index);

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

    bool isInside(const Vect &point) const {
        return point.x <= center.x + dimension.x / 2 &&
               point.x >= center.x - dimension.x / 2 &&
               point.y <= center.y + dimension.y / 2 &&
               point.y >= center.y - dimension.y / 2;
    }
};
