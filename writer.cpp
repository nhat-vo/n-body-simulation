#include "writer.hpp"
#include "common.hpp"

using namespace config;

void Writer::write(const std::vector<Vect> &r) {
    for (const Vect &v : r) {
        file << v.x - canvas_width/2 << " " << v.y - canvas_height/2<< " ";
    }
    file << std::endl;
}