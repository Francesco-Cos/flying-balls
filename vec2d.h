#ifndef VEC2D_H_INCLUDED
#define VEC2D_H_INCLUDED

#include <cmath>
#include <iostream>

class vec2d {
public:
    double x;
    double y;

    vec2d & operator = (const vec2d & other) {
	x = other.x;
	y = other.y;
	return *this;
    }

    vec2d & operator += (const vec2d & other) {
	x += other.x;
	y += other.y;
	return *this;
    }

    vec2d & operator -= (const vec2d & other) {
	x -= other.x;
	y -= other.y;
	return *this;
    }

    vec2d & operator *= (double l) {
	x *= l;
	y *= l;
	return *this;
    }

    vec2d operator + (const vec2d & other) const {
	return vec2d{x + other.x, y + other.y};
    }

    vec2d operator - (const vec2d & other) const {
	return vec2d{x - other.x, y - other.y};
    }

    vec2d operator * (double l) const {
	return vec2d{x*l, y*l};
    }

    vec2d operator / (double a) const {
	return vec2d{x/a, y/a};
    }

    static vec2d norm(const vec2d & v) {
        return v / module(v);
    } 

    static bool less(const vec2d & v1, const vec2d & v2) {
    if (v1.x >= 0 && v2.x < 0)
        return true;
    if (v1.x < 0 && v2.x >= 0)
        return false;
    if (v1.x == 0 && v2.x == 0) {
        if (v1.y >= 0 || v2.y >= 0)
            return v1.y > v2.y;
        return v2.y > v1.y;
    }

    // compute the cross product of vectors (center -> a) x (center -> b)
    int det = (v1.x) * (v2.y) - (v2.x) * (v1.y);
    if (det < 0)
        return true;
    if (det > 0)
        return false;

    // points a and b are on the same line from the center
    // check which point is closer to the center
    int d1 = (v1.x) * (v1.x) + (v1.y) * (v1.y);
    int d2 = (v2.x) * (v2.x) + (v2.y) * (v2.y);
    return d1 > d2;
    }

    vec2d & rotate (double angle) {
	double sin_a = sin(angle);
	double cos_a = cos(angle);
	double x1 = x;
	double y1 = y;
	x = cos_a*x1 - sin_a*y1;
	y = sin_a*x1 + cos_a*y1;
	return *this;
    }

    static double dot (const vec2d & a, const vec2d & b) {
	return a.x*b.x + a.y*b.y;
    }

    static double cross (const vec2d & a, const vec2d & b) {
	return a.x*b.y - a.y*b.x;
    }

    static double module (const vec2d & a) {
        return std::sqrt(dot(a, a));
    }

    friend std::ostream& operator <<(std::ostream&os, const vec2d& v){
    return os << "(" << v.x << ", " << v.y << ")";
}

};

static inline vec2d operator * (double l, const vec2d & v) {
    return vec2d{v.x*l, v.y*l};
}


#endif
