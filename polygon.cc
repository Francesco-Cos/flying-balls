#include <vector>
#include <cmath>
#include <cassert>
#include <iostream>
#include <random>

#include "polygon.h"
#include "vec2d.h"

polygon inner;
polygon outer;
unsigned int polygon_size = 5;
unsigned int track = 0;

void polygons_init()
{
    std::vector<vec2d> vec;
    for (int i = 0; i < polygon_size ; i++) 
    {
        vec2d temp{double(rand() % (width/2) + (-width/4)), double(rand() % (height/2) + (-height/4))};
        if (vec2d::module(temp) < 100) {
            vec2d norm = vec2d::norm(temp);
            temp = norm * 150; 
        }
        vec.push_back(temp);
    }
    
    inner = polygon{vec};
    inner.set_center();
    std::sort(inner.points.begin(), inner.points.end(), vec2d::less);

    outer = create_outer(inner);
    outer.set_center();
    std::sort(outer.points.begin(), outer.points.end(), vec2d::less);

    for (auto& point : inner.points) {
        point.x += width/2;
        point.y += height/2;
    }

    for (auto& point : outer.points) {
        point.x += width/2;
        point.y += height/2;
    }
    for (auto& point : inner.points) {
        std::cout << "point: " << point << std::endl; 
    }
}

polygon create_outer(polygon& p)
{
    std::vector<vec2d> points;
    points.reserve(p.points.size());
    for (int i = 1; i <= p.points.size(); i++)
    {
        vec2d& point = p.points[i%p.points.size()];
        vec2d vec1 = vec2d::norm(point - p.points[(i-1)%p.points.size()]);
        vec2d vec2 = vec2d::norm(point - p.points[(i+1)%p.points.size()]);
        vec2d dir = (vec1 + vec2);
        if (vec2d::module(dir) < 1) dir = vec2d::norm(dir);
        if (vec2d::dot(dir, point) < 0) dir *= -1;

        points.push_back(point + dir * 80);
    }
        polygon ret{points};

    return ret;
}

void polygon::draw(cairo_t *cr) const
{
    cairo_save(cr);
    cairo_new_path(cr);
    for (auto &point : points)
        cairo_line_to(cr, point.x, point.y);
    cairo_line_to(cr, points[0].x, points[0].y);
    cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
    cairo_set_line_width(cr, 1.0);
    cairo_stroke(cr);
    cairo_restore(cr);
}

void polygon_draw(cairo_t *cr)
{
    inner.draw(cr);
    outer.draw(cr);
}