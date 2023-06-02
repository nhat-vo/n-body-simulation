#pragma once

#include "common.hpp"
#include "visualizer.hpp"

using namespace config;

#ifdef VISUALIZE

void single_thread(Scenario &bodies, size_t n_threads, Drawer &drawer);
void multi_thread_1(Scenario &bodies, size_t n_threads, Drawer &drawer);
void multi_thread_2(Scenario &bodies, size_t n_threads, Drawer &drawer);
void barnes_hut(Scenario &bodies, size_t n_threads, Drawer &drawer);

#else

void single_thread(Scenario &bodies, size_t n_threads);
void multi_thread_1(Scenario &bodies, size_t n_threads);
void multi_thread_2(Scenario &bodies, size_t n_threads);
void barnes_hut(Scenario &bodies, size_t n_threads);

#endif