#include "writer.hpp"

void Writer::write(const std::vector<Vect> &r) {
    for (const Vect &v : r) {
        file << v.x << " " << v.y << " ";
    }
    file << std::endl;
}