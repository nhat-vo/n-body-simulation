#pragma once

#include "common.hpp"
#include "visualizer.hpp"
#include "writer.hpp"

using namespace config;

void initialize_bodies(Scenario &bodies, size_t n_bodies, std::vector<std::string> &colors);

#ifdef VISUALIZE

void single_thread(Scenario &bodies, size_t n_threads, Drawer &drawer);
void multi_thread_1(Scenario &bodies, size_t n_threads, Drawer &drawer);
void multi_thread_2(Scenario &bodies, size_t n_threads, Drawer &drawer);
void barnes_hut(Scenario &bodies, size_t n_threads, Drawer &drawer);
void barnes_hut_multi(Scenario &bodies, size_t n_threads, Drawer &drawer);

#elif defined(WRITE)

void single_thread(Scenario &bodies, size_t n_threads, Writer &writer);
void multi_thread_1(Scenario &bodies, size_t n_threads, Writer &writer);
void multi_thread_2(Scenario &bodies, size_t n_threads, Writer &writer);
void barnes_hut(Scenario &bodies, size_t n_threads, Writer &writer);
void barnes_hut_multi(Scenario &bodies, size_t n_threads, Writer &writer);

#else

void single_thread(Scenario &bodies, size_t n_threads);
void multi_thread_1(Scenario &bodies, size_t n_threads);
void multi_thread_2(Scenario &bodies, size_t n_threads);
void barnes_hut(Scenario &bodies, size_t n_threads);
void barnes_hut_multi(Scenario &bodies, size_t n_threads);

#endif
