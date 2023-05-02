#ifndef BALLS_H_INCLUDED
#define BALLS_H_INCLUDED

#include <gtk/gtk.h>

#include "vec2d.h"

class ball_face;

class ball {
public:
    unsigned int border;
    unsigned int radius;
    vec2d position;
    vec2d velocity;
    vec2d force;

    double angle;
    double v_angle;

    ball_face * face;

    void draw (cairo_t * cr) const;
};

extern ball * balls;
extern unsigned int n_balls;
extern unsigned int border_particles;
extern unsigned int fluid;
extern unsigned int show_fluid;
extern unsigned int rc;
extern unsigned int max_rep;
extern double sigma;

extern unsigned int radius_min;
extern unsigned int radius_max;
extern unsigned int radius_particle;

extern unsigned int v_max;
extern unsigned int v_min;
extern double eps;


extern unsigned int v_angle_min;
extern unsigned int v_angle_max;

extern const char * face_filename;
extern int face_rotation;

extern void balls_init ();
extern void balls_destroy ();
extern void ball_update_pos (ball * p);
extern void ball_update_state (ball * p);
extern void ball_update_state_fluid (ball * p);
extern vec2d ball_calculate_force (ball * p);
extern vec2d ball_calculate_distance (vec2d d1, vec2d d2);
extern void ball_ball_collision (ball * p, ball * q);
extern void ball_reposition (ball * b);
extern void balls_draw (cairo_t * cr);
extern void shoot_draw (cairo_t * cr);

extern void restitution_coefficient_draw (cairo_t * cr);
extern void restitution_coefficient_set (double c);
extern double restitution_coefficient_get ();
extern void restitution_coefficient_change (double d);

#endif
