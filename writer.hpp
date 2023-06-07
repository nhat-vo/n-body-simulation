#pragma once

#ifndef WRITER
#define WRITER

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "vect.hpp"

class Writer {
    private:
        std::ofstream file;
        std::string filename;

    public:
        Writer(const std::string &filename) : filename(filename) {
            file.open(filename);
            if (!file.is_open()) {
                std::cout << "Cannot open file " << filename << std::endl;
                exit(1);
            }
        }

        ~Writer() { file.close(); }

        void write(const std::vector<Vect> &r);
};

#endif