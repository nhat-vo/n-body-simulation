#include "common.hpp"
#include "vect.hpp"
#include "visualizer.hpp"

#include <chrono>
#include <cmath>
#include <iostream>
#include <list>
#include <thread>
#include <vector>
#include <tuple>

using namespace config;

struct Body {
    double m;
    Vect r;
    Vect v;
    std::string color;
    Body(double m, Vect r, Vect v, std::string color) : m(m), r(r), v(v), color(color) {}
};

struct QuadNode {
    QuadNode *children[4]{nullptr, nullptr, nullptr, nullptr};
    double mass;
    Vect center_of_mass;
    Vect TopLeft;
    Vect dimension;
    std::vector<Body *> bodies;
    QuadNode(Vect TopLeft, Vect dimension) : TopLeft(TopLeft), dimension(dimension) {
        bodies = std::vector<Body *>();
        Vect center_of_mass = Vect(TopLeft.x + dimension.x / 2, TopLeft.y + dimension.y / 2);
        double mass = 0;
    }
};

/* ---------------- BARNES-HUT AUX ------------------- */

void updateCenterOfMass(QuadNode *curr, Body *body) {
    double m1, m2, m;

    m1 = curr->mass;
    m2 = body->m;
    m = m1 + m2;

    curr->center_of_mass = (curr->center_of_mass * m1 + body->r * m2) / m;
    curr->mass = m;
}

int getQuad(Vect TopLeft, Vect dimension, Vect r) {
    Vect center = TopLeft + Vect(dimension.x / 2, - dimension.y / 2);
    return r.x < center.x ? (r.y < center.y ? 2 : 0)
                            : (r.y < center.y ? 3 : 1);
}

Vect calcTopLeft(Vect TopLeft, Vect dimension, int quad) {
    Vect center = TopLeft + Vect(dimension.x / 2, - dimension.y / 2);
    switch (quad) {
        case 0: return Vect(TopLeft.x, TopLeft.y);
        case 1: return Vect(center.x, TopLeft.y);
        case 2: return Vect(TopLeft.x, center.y);
        case 3: return Vect(center.x, center.y);
        default: return Vect(0, 0);
    };
}

/* ---------------- BARNES-HUTT ------------------- */

void BarnesHutRecursion(QuadNode *curr, Body *body) {

    size_t n = curr->bodies.size();
    curr->bodies.push_back(body);

    // If node is empty, add body to node
    if (n == 0) {
        updateCenterOfMass(curr, body);
        return;
    }

    /* If there is only one body in the node,
     * create a new define the sub node of the old
     * because the node was not split yet
     */
    if (n == 1) {
        int quad_old_body = getQuad(curr->TopLeft, curr->dimension, curr->bodies[0]->r);
        Vect QuadTopLeft = calcTopLeft(curr->TopLeft, curr->dimension, quad_old_body);
        curr->children[quad_old_body] = new QuadNode(QuadTopLeft, curr->dimension / 2);
        BarnesHutRecursion(curr->children[quad_old_body], curr->bodies[0]);
    }

    /* Get quad where new body belongs */
    int quad_new_body = getQuad(curr->TopLeft, curr->dimension, body->r);

    /* If the new node is empty, initialize the node */
    if (curr->children[quad_new_body] == nullptr) {
        Vect QuadTopLeft = calcTopLeft(curr->TopLeft, curr->dimension, quad_new_body);
        curr->children[quad_new_body] = new QuadNode(QuadTopLeft, curr->dimension / 2);
    }
    /* Add the body to the node and update center of mass */
    BarnesHutRecursion(curr->children[quad_new_body], body);
    updateCenterOfMass(curr, body);
}

bool isFarEnough(QuadNode *curr, Body *body) {
    Vect dr = curr->center_of_mass - body->r;
    double dist_sq = std::max(dr.norm2(), 1e-6);
    double dist = sqrt(dist_sq);
    return curr->dimension.x / dist < theta;
}

void constructBarnesHutTree(Scenario bodies) {

    // Construct root with first body
    QuadNode *root = new QuadNode(Vect(0, canvas_height), Vect(canvas_width, canvas_height));
    root->mass = bodies.m[0];
    root->center_of_mass = bodies.r[0];
    Body b = Body(bodies.m[0], bodies.r[0], bodies.v[0], bodies.colors[0]);
    root->bodies.push_back(&b);

    // Construct tree with remaining bodies
    for (int i = 1; i < bodies.r.size(); i++) {
        Body *b = new Body(bodies.m[i], bodies.r[i], bodies.v[i], bodies.colors[i]);
        BarnesHutRecursion(root, b);
    }

    /* Calculate the force exerted */
    for (int i = 0; i < bodies.r.size(); i++) {
        Body b = Body(bodies.m[i], bodies.r[i], bodies.v[i], bodies.colors[i]);
        std::list<QuadNode *> stack;
        stack.push_back(root);
        while (!stack.empty()) {
            QuadNode *curr = stack.back();
            stack.pop_back();
            size_t n = curr->bodies.size();
            if (n == 1) {
                if ((curr->bodies[0]->r - b.r).norm2() > 0) {
                    Vect dr = curr->bodies[0]->r - bodies.r[i];
                    double dist_sq = std::max(dr.norm2(), 1e-6);
                    double dist = sqrt(dist_sq);
                    double force_mag = G * curr->bodies[0]->m * b.m / dist_sq;
                    Vect force = dr * (force_mag / dist);
                    bodies.v[i] += force * (dt / b.m);
                }
            } else if (isFarEnough(curr, &b)) {
                Vect dr = curr->center_of_mass - bodies.r[i];
                double dist_sq = std::max(dr.norm2(), 1e-6);
                double dist = sqrt(dist_sq);
                double force_mag = G * curr->mass * b.m / dist_sq;
                Vect force = dr * (force_mag / dist);
                bodies.v[i] += force * (dt / b.m);
            } else {
                for (int j = 0; j < 4; j++) {
                    if (curr->children[j])
                        stack.push_back(curr->children[j]);
                }
            }
        }
    }
}