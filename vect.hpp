#pragma once
#include <math.h>
#include <iostream>

class Vect {
  public:
    double x, y;
    Vect(double x, double y) : x(x), y(y) {}
    Vect() : x(0), y(0) {}

    double norm2() { return x * x + y * y; }
    double norm() { return sqrt(x * x + y * y); }

    void operator+=(const Vect &b) {
        x += b.x;
        y += b.y;
    }
    void operator-=(const Vect &b) {
        x -= b.x;
        y -= b.y;
    }
    Vect operator+(const Vect &b) const { return {x + b.x, y + b.y}; }
    Vect operator-(const Vect &b) const { return {x - b.x, y - b.y}; }
    Vect operator*(const Vect &b) const { return {x * b.x, y * b.y}; }
    Vect operator/(const Vect &b) const { return {x / b.x, y / b.y}; }

    Vect operator+(double b) { return {x + b, y + b}; }
    Vect operator-(double b) { return {x - b, y - b}; }
    Vect operator*(double b) { return {x * b, y * b}; }
    Vect operator/(double b) { return {x / b, y / b}; }
};
std::ostream &operator<<(std::ostream &os, const Vect &v) {
    os << "(" << v.x << ", " << v.y << ")";
    return os;
}