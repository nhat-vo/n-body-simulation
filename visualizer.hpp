/*
 * Licensed under the ImageMagick License (the "License"); you may not use
 * this file except in compliance with the License.  You may obtain a copy of
 * the License at https://imagemagick.org/script/license.php Unless required by
 * applicable law or agreed to in writing, software distributed under the
 * License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS
 * OF ANY KIND, either express or implied.  See the License for the specific
 * language governing permissions and limitations under the License.
 *
 *
 *
 * This is a simple Drawer class using Magick++ for visualization. In order to
 * compile with this, you need to have Magick++ installed, which comes with the
 * imagemagick package. https://imagemagick.org/. You need to add -DVISUALIZE
 * flag in compilation. Have a look at the Makefile for more information. If you
 * cannot install this, simply exclude the flag.
 *
 * This Drawer class is designed in a call-and-forget manner to make it simple
 * for integration. You just need to create a Drawer object, and call
 * trigger_draw() ONCE per timestep. The image will be automatically exported on
 * destruction. Have a look at the main file to see how this is used.
 *
 */
#pragma once

#ifdef VISUALIZE

#include "common.hpp"
#include "vect.hpp"
#include <Magick++.h>
#include <iostream>
#include <thread>
using namespace config;

class Drawer {
  private:
    const Magick::Geometry image_size;
    const Magick::Color background_color;
    const std::vector<std::string> &colors;

    std::list<Magick::Image> images;
    double draw_t = 0;
    std::vector<Vect> draw_r;

    std::thread draw_thread;

    static void draw(Drawer *drawer_pt);

  public:
    Drawer(const std::vector<std::string> &colors);
    ~Drawer();
    void trigger_draw(double t, std::vector<Vect> *curr_r);
};

#endif
