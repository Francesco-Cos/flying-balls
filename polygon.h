#ifndef POLYGON_H_INCLUDED
#define POLYGON_H_INCLUDED

#include "balls.h"

#include "vec2d.h"
#include "game.h"
#include <vector>

class polygon {
public:
    std::vector<vec2d> points;

    void set_center() {
        vec2d centroid{.x = 0, .y = 0};
        for(auto& point : points) 
            centroid += point;
        centroid = centroid / points.size();
        std::cout << "centroid: " << centroid << std::endl;
        for(auto& point : points) {
            point.x += (- centroid.x);
            point.y += (- centroid.y);
        }
    }

    void draw (cairo_t * cr) const;
};

extern unsigned int track;

extern polygon inner;
extern polygon outer;
extern unsigned int polygon_size;
extern void polygons_init();
extern polygon create_outer(polygon& p);
extern void polygon_draw (cairo_t * cr);

#endif