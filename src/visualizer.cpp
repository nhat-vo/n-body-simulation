#include "../headers/visualizer.hpp"

#ifdef VISUALIZE

Drawer::Drawer(const std::vector<std::string> &colors)
    : image_size(canvas_width, canvas_height), background_color("white"),
      colors(colors) {}

Drawer::~Drawer() {
    if (draw_thread.joinable())
        draw_thread.join();
    std::cout << "Exporting images to image.gif..." << std::endl;
    Magick::writeImages(images.begin(), images.end(), "image.gif");
    std::cout << "Image exported to image.gif" << std::endl;
}

void Drawer::trigger_draw(double t, std::vector<Vect> *curr_r) {
    if (t >= draw_t) {
        // wait until the draw thread has finished drawing
        if (draw_thread.joinable())
            draw_thread.join();

        draw_t += draw_dt;
        draw_r = *curr_r;

        draw_thread = std::thread(Drawer::draw, this);
    }
}

void Drawer::draw(Drawer *drawer_pt) {
    // auto start = std::chrono::steady_clock::now();
    Drawer &drawer = *drawer_pt;

    drawer.images.emplace_back(drawer.image_size, drawer.background_color);
    auto &frame = drawer.images.back();

    for (size_t i = 0; i < drawer.draw_r.size(); ++i) {
        const Vect &r = drawer.draw_r[i];
        double x = r.x, y = (double)canvas_height - r.y;

        frame.fillColor(drawer.colors[i]);
        frame.draw(Magick::DrawableCircle(x, y, x - 5, y));
    }

    frame.strokeColor("black");
    frame.draw(Magick::DrawableText(
        10, canvas_height - 10,
        "t = " + std::to_string((int)(drawer.draw_t - draw_dt))));
    // auto end = std::chrono::steady_clock::now();
    // std::chrono::duration<double> duration = end - start;
    // std::cout << "draw time: " << duration.count() << std::endl;
}

#endif
