#include <vector>
#include <cmath>
#include <cassert>
#include <iostream>
#include <random>

#include "polygon.h"
#include "game.h"
#include "balls.h"
#include "gravity.h"
#include "cell_list.h"

unsigned int radius_min = 5;
unsigned int radius_max = 10;
unsigned int radius_particle = 6;
unsigned int max_rep = 35;
unsigned int rc = radius_particle * 2;
unsigned int border_velocity = 100;
double sigma = 1.0;

unsigned int v_max = 100;
unsigned int v_min = 0;
double eps = 1.0;

unsigned int v_angle_min = 0;
unsigned int v_angle_max = 100;

vec2d o{.x = 0, .y = 0};
std::default_random_engine generator;
std::normal_distribution<double> distribution(0, 1.0);

ball *balls = nullptr;
unsigned int n_balls = 50;
unsigned int border_particles = 5;
unsigned int fluid = 0;
unsigned int show_fluid = 1;

static double C_r = 1.0;

void restitution_coefficient_set(double c)
{
	C_r = c;
	if (C_r > 1.0)
		C_r = 1.0;
	else if (C_r < 0.0)
		C_r = 0.0;
}

double restitution_coefficient_get()
{
	return C_r;
}

void restitution_coefficient_change(double d)
{
	C_r += d;
	if (C_r > 1.0)
		C_r = 1.0;
	else if (C_r < 0.0)
		C_r = 0.0;
}

void restitution_coefficient_draw(cairo_t *cr)
{
	static const double margin = 20;
	cairo_save(cr);
	cairo_new_path(cr);
	cairo_move_to(cr, margin, margin);
	cairo_line_to(cr, margin + (width - 2 * margin) * C_r, margin);
	cairo_set_source_rgb(cr, 1.0, 1.0, 0.0);
	cairo_set_line_width(cr, margin / 2);
	cairo_stroke(cr);
	cairo_restore(cr);
}

static vec2d random_velocity()
{
	double r2;
	vec2d v;
	do
	{
		v.x = v_min + rand() % (v_max + 1 - v_min);
		v.y = v_min + rand() % (v_max + 1 - v_min);
		r2 = vec2d::dot(v, v);
	} while (r2 > v_max * v_max || r2 < v_min * v_min);
	if (rand() % 2)
		v.x = -v.x;
	if (rand() % 2)
		v.y = -v.y;
	return v;
}

void balls_init_state()
{
	static const int border = 10;
	int w = width < 2 * border ? 1 : width - 2 * border;
	int h = height < 2 * border ? 1 : height - 2 * border;

	for (unsigned int i = 0; i < n_balls - border_particles; ++i)
	{
		balls[i].position.x = border + rand() % w;
		balls[i].position.y = border + rand() % h;
		balls[i].velocity = fluid ? o : random_velocity();
		balls[i].border = 0;
		balls[i].segment = -1;
		balls[i].inner = -1;
		balls[i].radius = fluid ? radius_particle : radius_min + rand() % (radius_max + 1 - radius_min);
		unsigned int v_angle_360 = (v_angle_min + rand() % (v_angle_max + 1 - v_angle_min)) % 360;
		balls[i].v_angle = 2 * M_PI * v_angle_360 / 360;
		balls[i].angle = (rand() % 360) * 2 * M_PI / 360;
	}
	if (fluid)
	{
		for (unsigned int i = n_balls - border_particles; i < n_balls; ++i)
		{
			balls[i].position.x = i % 2 ? border : width - border;
			balls[i].position.y = border + rand() % h;
			balls[i].velocity.x = 0;
			balls[i].velocity.y = i % 2 ? 300 : -300;
			balls[i].border = 1;
			balls[i].segment = -1;
			balls[i].inner = -1;
			balls[i].radius = radius_particle;
			balls[i].v_angle = 0;
			balls[i].angle = 0;
		}

		for (unsigned int i = 0; i < n_balls - border_particles; ++i)
		{
			balls[i].force = ball_calculate_force(&balls[i]);
		}
	}
}

std::vector<int> weights() {
	int particles = n_balls - 2 * inner.points.size();
	std::vector<double> segments = {};
	for (int i = 0; i < inner.points.size(); i++) {
		vec2d vec_inner = inner.points[i];
		int next_vec = i + 1 >= inner.points.size() ? 0 : i + 1;
		vec2d next_vec_inner = inner.points[next_vec];
		vec2d seg = next_vec_inner - vec_inner;
		double seg_len = vec2d::module(seg);
		segments.push_back(seg_len);
	}
	double perimeter = 0;
	for (auto seg : segments) perimeter += seg;
	std::vector<int> weights = {};
	for (int i = 0; i < inner.points.size(); i++) {
		weights.push_back(floor((segments[i]/perimeter)*particles));
		std::cout << weights[i] << std::endl;
	}
	return weights;
}

void balls_init_state_track()
{
	std::vector<int> w = weights();
	std::cout << w.size() << std::endl;
	int start = 0;
	for (int j = 0; j < inner.points.size(); j++) {
		std::cout << inner.points[j] << std::endl;
	for (unsigned int i = start; i < start + w[j]; ++i)
	{
		vec2d vec_inner = inner.points[j];
		int next_vec = j + 1 >= inner.points.size() ? 0 : j + 1;
		vec2d next_vec_inner = inner.points[next_vec];
		vec2d seg = next_vec_inner - vec_inner;
		vec2d vec_outer = outer.points[j];
		vec2d t = vec_outer - vec_inner;
		balls[i].position = vec_inner + (t * (30 + rand() % 40) / 100) + (seg * (rand() % 100) / 100);
		balls[i].velocity = fluid ? o : random_velocity();
		balls[i].border = 0;
		balls[i].segment = -1;
		balls[i].inner = -1;
		balls[i].radius = fluid ? radius_particle : radius_min + rand() % (radius_max + 1 - radius_min);
		unsigned int v_angle_360 = (v_angle_min + rand() % (v_angle_max + 1 - v_angle_min)) % 360;
		balls[i].v_angle = 2 * M_PI * v_angle_360 / 360;
		balls[i].angle = (rand() % 360) * 2 * M_PI / 360;
	}
	start += w[j];
	}

	for (unsigned int i = n_balls - 2 * inner.points.size(); i < n_balls - inner.points.size(); ++i)
	{
		int point = i % (n_balls - 2 * inner.points.size());
		balls[i].position = inner.points[point];
		int next_point = point + 1 >= inner.points.size() ? 0 : point + 1;
		vec2d vel_dir = vec2d::norm(inner.points[next_point] - inner.points[point]);
		balls[i].velocity = vel_dir * border_velocity;
		balls[i].border = 1;
		balls[i].segment = point;
		balls[i].inner = 1;
		balls[i].radius = fluid ? radius_particle : radius_min + rand() % (radius_max + 1 - radius_min);
		unsigned int v_angle_360 = (v_angle_min + rand() % (v_angle_max + 1 - v_angle_min)) % 360;
		balls[i].v_angle = 2 * M_PI * v_angle_360 / 360;
		balls[i].angle = (rand() % 360) * 2 * M_PI / 360;
	}

	for (unsigned int i = n_balls - inner.points.size(); i < n_balls; ++i)
	{
		int point = i % (n_balls - inner.points.size());
		balls[i].position = outer.points[point];
		int next_point = point + 1 >= inner.points.size() ? 0 : point + 1;
		vec2d vel_dir = vec2d::norm(outer.points[next_point] - outer.points[point]);
		balls[i].velocity = vel_dir * border_velocity;
		balls[i].border = 1;
		balls[i].segment = point;
		balls[i].inner = 0;
		balls[i].radius = fluid ? radius_particle : radius_min + rand() % (radius_max + 1 - radius_min);
		unsigned int v_angle_360 = (v_angle_min + rand() % (v_angle_max + 1 - v_angle_min)) % 360;
		balls[i].v_angle = 2 * M_PI * v_angle_360 / 360;
		balls[i].angle = (rand() % 360) * 2 * M_PI / 360;
	}
}

void ball_top_bottom_collision(ball *p)
{
	if (p->position.y + p->radius > height)
	{ /* bottom wall */
		if (p->velocity.y > 0)
		{
			p->position.y -= p->position.y + p->radius - height;
			p->velocity.y = -C_r * p->velocity.y;
		}
	}
	else if (p->position.y < p->radius)
	{ /* top wall */
		if (p->velocity.y < 0)
		{
			p->position.y += p->radius - p->position.y;
			p->velocity.y = -C_r * p->velocity.y;
		}
	}
}

void ball_top_bottom_collision_fluid(ball *p)
{
	if (p->position.y > height)
	{ /* bottom wall */
		if (p->velocity.y > 0)
		{
			p->position.y = 0;
		}
	}
	else if (p->position.y < 0)
	{ /* top wall */
		if (p->velocity.y < 0)
		{
			p->position.y = height;
		}
	}
}

void ball_polygon_collision(ball *b, polygon *p)
{
	int size = p->points.size();
	for (int i = 1; i <= size; i++)
	{
		vec2d &point = p->points[i % size];
		vec2d point2 = p->points[(i + 1) % size];
		vec2d vec1 = point2 - point;
		vec2d vec2 = b->position - point;
		vec2d vec3 = b->position - point2;
		vec2d projection = vec1 * (vec2d::dot(vec1, vec2) / vec2d::dot(vec1, vec1));
		vec2d normal = vec2d::norm(vec2 - projection);
		vec2d reflected_velocity = b->velocity - 2 * vec2d::dot(b->velocity, normal) * normal;
		double dist = vec2d::module(vec2 - projection);
		if (dist < radius_particle && vec2d::module(vec2) < vec2d::module(vec1) && vec2d::module(vec3) < vec2d::module(vec1))
		{
			b->position += normal;
			b->velocity = reflected_velocity;
		}
	}
}

void ball_walls_collision(ball *p)
{
	if (p->position.x + p->radius > width)
	{ /* right wall */
		if (p->velocity.x > 0)
		{
			p->position.x -= p->position.x + p->radius - width;
			p->velocity.x = -C_r * p->velocity.x;
		}
	}
	else if (p->position.x < p->radius)
	{ /* left wall */
		if (p->velocity.x < 0)
		{
			p->position.x += p->radius - p->position.x;
			p->velocity.x = -C_r * p->velocity.x;
		}
	}
	if (fluid)
	{
		ball_top_bottom_collision_fluid(p);
	}
	else
	{
		ball_top_bottom_collision(p);
	}
}

vec2d ball_calculate_force(ball *p)
{
	vec2d Fc{.x = 0, .y = 0};
	vec2d Fd{.x = 0, .y = 0};
	vec2d Fr{.x = 0, .y = 0};
	if (!clist)
	{
		// std::cout << "wrong" << std::endl;
		for (unsigned int i = 0; i < n_balls; ++i) {
			double n = distribution(generator);
			vec2d sub_p = p->position - balls[i].position;
			vec2d sub_v = p->velocity - balls[i].velocity;
			double rij = vec2d::module(sub_p);
			double wr = rc - rij;
			double wd = wr * wr;
			vec2d Rij = rij < eps ? o : sub_p / rij;
			Fc += wr > 0 ? max_rep * wr * Rij : o;
			Fd += wr > 0 ? -(sigma * sigma) / 2 * wd * vec2d::dot(Rij, sub_v) * Rij : o;
			Fr += wr > 0 ? sigma * wr * n * Rij : o;
		}
	}
	else
	{
		std::vector<ball*> nh = neighbourhood(p);
		for (ball * ball : nh)
		{
			double n = distribution(generator);
			vec2d sub_p = p->position - ball->position;
			vec2d sub_v = p->velocity - ball->velocity;
			double rij = vec2d::module(sub_p);
			double wr = rc - rij;
			double wd = wr * wr;
			vec2d Rij = rij < eps ? o : sub_p / rij;
			Fc += wr > 0 ? max_rep * wr * Rij : o;
			Fd += wr > 0 ? -(sigma * sigma) / 2 * wd * vec2d::dot(Rij, sub_v) * Rij : o;
			Fr += wr > 0 ? sigma * wr * n * Rij : o;
		}
	}
	return isnan(Fc.x) ? o : Fc + Fd + Fr;
}

void ball_update_pos(ball *p)
{
	if (fluid)
	{
		if (!p->border)
		{
			p->position += delta * p->velocity + delta * delta * p->force / 2.0;
		}
	}
}

void ball_update_state(ball *p)
{
	vec2d g = gravity_vector(p);
	p->position += delta * p->velocity + delta * delta * g / 2.0;
	p->velocity += delta * g;
	p->angle += delta * p->v_angle;
	while (p->angle >= 2 * M_PI)
		p->angle -= 2 * M_PI;
	while (p->angle < 0)
		p->angle += 2 * M_PI;
	ball_walls_collision(p);
}

void ball_update_state_fluid(ball *p)
{
	if (!p->border)
	{
		vec2d old_v = p->velocity;
		p->velocity += delta * p->force / 2;
		vec2d old_f = p->force;
		p->force = ball_calculate_force(p);
		p->velocity = old_v + delta * (p->force + old_f) / 2;
		if (vec2d::module(p->velocity) > 100) p->velocity = vec2d::norm(p->velocity) * 100;
		p->angle += delta * p->v_angle;
		while (p->angle >= 2 * M_PI)
			p->angle -= 2 * M_PI;
		while (p->angle < 0)
			p->angle += 2 * M_PI;
	}
	else
	{
		p->position += delta * p->velocity;
		int next_segment = p->segment + 1 >= inner.points.size() ? 0 : p->segment + 1;
		if (track && p->inner &&
			p->position.x <= inner.points[next_segment].x + 1.5 &&
			p->position.x >= inner.points[next_segment].x - 1.5 &&
			p->position.y <= inner.points[next_segment].y + 1.5 &&
			p->position.y >= inner.points[next_segment].y - 1.5)
		{
			p->segment = next_segment;
			p->position = inner.points[next_segment];
			p->velocity = vec2d::norm(inner.points[(next_segment + 1) % inner.points.size()] - inner.points[next_segment]) * border_velocity;
		}
		else if (track && !p->inner &&
				 p->position.x <= outer.points[next_segment].x + 1.5 &&
				 p->position.x >= outer.points[next_segment].x - 1.5 &&
				 p->position.y <= outer.points[next_segment].y + 1.5 &&
				 p->position.y >= outer.points[next_segment].y - 1.5)
		{
			p->segment = next_segment;
			p->position = outer.points[next_segment];
			p->velocity = vec2d::norm(outer.points[(next_segment + 1) % inner.points.size()] - outer.points[next_segment]) * border_velocity;
		}
		if (track)
			p->velocity = vec2d::norm(p->velocity) * border_velocity;
	}
	if (track && !p->border)
	{
		ball_polygon_collision(p, &inner);
		ball_polygon_collision(p, &outer);
	}
	else if (!track)
	{
		ball_walls_collision(p);
	}
}

void ball_ball_collision(ball *p, ball *q)
{
	if (!fluid)
	{
		vec2d pq = q->position - p->position;
		double d2 = vec2d::dot(pq, pq);
		double r = p->radius + q->radius;
		if (d2 <= r * r)
		{
			vec2d pq_v = q->velocity - p->velocity;

			double mp = p->radius * p->radius;
			double mq = q->radius * q->radius;
			double m_total = mp + mp;

			double d = sqrt(d2);
			vec2d pq_overlap = (r - d) / d * pq;
			p->position -= pq_overlap * mq / m_total;
			q->position += q->border ? o : pq_overlap * mp / m_total;

			double f = vec2d::dot(pq_v, pq);

			if (f < 0)
			{
				f /= d2 * (mp + mq);
				p->velocity += 2 * C_r * mq * f * pq;
				q->velocity -= q->border ? o : 2 * C_r * mp * f * pq;
			}
		}
	}
}

void ball_reposition(ball *b)
{
	if (b->position.x < b->radius)
		b->position.x = b->radius;
	else if (b->position.x + b->radius > width)
		b->position.x = width - b->radius;
	if (b->position.y < b->radius)
		b->position.y = b->radius;
	else if (b->position.y + b->radius > height)
		b->position.y = height - b->radius;
}

const char *face_filename = 0;
int face_rotation = 0;

static const double linear_rotation_unit = 2.0;

static std::vector<ball_face *> faces;

class ball_face
{
public:
	ball_face(unsigned int radius, cairo_surface_t *face, int rotation, ball *b);
	~ball_face();
	cairo_surface_t *get_surface(double angle) const;

private:
	unsigned int rotations;
	std::vector<cairo_surface_t *> c_faces;
};

static double random_color_component()
{
	return 1.0 * (rand() % 200 + 56) / 255;
};

ball_face::ball_face(unsigned int radius, cairo_surface_t *face, int rotation, ball *b)
{
	if (face && rotation)
	{
		rotations = 2 * M_PI * radius / linear_rotation_unit;
	}
	else
	{
		rotations = 1;
	}
	c_faces.resize(rotations);
	for (unsigned int i = 0; i < rotations; ++i)
	{
		c_faces[i] = gdk_window_create_similar_surface(gtk_widget_get_window(canvas),
													   CAIRO_CONTENT_COLOR_ALPHA,
													   2 * radius, 2 * radius);
		assert(c_faces[i]);
		cairo_t *ball_cr = cairo_create(c_faces[i]);
		cairo_translate(ball_cr, radius, radius);
		cairo_arc(ball_cr, 0.0, 0.0, radius, 0, 2 * M_PI);
		cairo_clip(ball_cr);

		if (face)
		{
			int face_x_offset = cairo_image_surface_get_width(face) / 2;
			int face_y_offset = cairo_image_surface_get_height(face) / 2;
			cairo_rotate(ball_cr, i * 2 * M_PI / rotations);
			cairo_scale(ball_cr, 1.0 * radius / face_x_offset, 1.0 * radius / face_y_offset);
			cairo_set_source_surface(ball_cr, face, -face_x_offset, -face_y_offset);
			cairo_paint(ball_cr);
		}
		else
		{
			cairo_pattern_t *pat;
			pat = cairo_pattern_create_radial(-0.2 * radius, -0.2 * radius, 0.2 * radius,
											  -0.2 * radius, -0.2 * radius, 1.3 * radius);

			double col_r = b->border ? 0.9 : 0.1; // random_color_component();
			double col_g = b->border ? 0.5 : 0.5; // random_color_component();
			double col_b = b->border ? 0.1 : 0.9; // random_color_component();
			cairo_pattern_add_color_stop_rgba(pat, 0, col_r, col_g, col_b, 1);
			cairo_pattern_add_color_stop_rgba(pat, 1, col_r / 3, col_g / 3, col_b / 3, 1);
			cairo_set_source(ball_cr, pat);
			cairo_arc(ball_cr, 0.0, 0.0, radius, 0, 2 * M_PI);
			cairo_fill(ball_cr);
			cairo_pattern_destroy(pat);
		}
		cairo_surface_flush(c_faces[i]);
		cairo_destroy(ball_cr);
	}
}

ball_face::~ball_face()
{
	for (auto f : c_faces)
		cairo_surface_destroy(f);
}

cairo_surface_t *ball_face::get_surface(double angle) const
{
	unsigned int face_id;
	if (rotations == 1)
		face_id = 0;
	else
	{
		face_id = rotations * angle / (2 * M_PI);
		assert(face_id < rotations);
		if (face_id >= rotations)
			face_id %= rotations;
	}
	return c_faces[face_id];
}

static void balls_init_faces()
{
	cairo_surface_t *face_surface = 0;

	if (face_filename)
	{
		face_surface = cairo_image_surface_create_from_png(face_filename);
		if (cairo_surface_status(face_surface) != CAIRO_STATUS_SUCCESS)
		{
			cairo_surface_destroy(face_surface);
			face_surface = 0;
			fprintf(stderr, "could not create surface from PNG file %s\n", face_filename);
		}
	}
	if (face_surface)
	{
		faces.assign(radius_max + 1 - radius_min, nullptr);
		for (ball *b = balls; b != balls + n_balls; ++b)
		{
			unsigned int r_idx = b->radius - radius_min;
			if (faces[r_idx] == nullptr)
				faces[r_idx] = new ball_face(b->radius, face_surface, face_rotation, b);
			b->face = faces[r_idx];
		}
		cairo_surface_destroy(face_surface);
	}
	else
	{
		faces.resize(n_balls);
		for (unsigned int i = 0; i < n_balls; ++i)
			balls[i].face = faces[i] = new ball_face(balls[i].radius, 0, face_rotation, &balls[i]);
	}
}

void ball::draw(cairo_t *cr) const
{
	if (show_fluid)
	{
		cairo_save(cr);
		cairo_translate(cr, position.x - radius, position.y - radius);
		cairo_set_source_surface(cr, face->get_surface(angle), 0, 0);
		cairo_paint(cr);
		cairo_restore(cr);
	}
}

void balls_draw(cairo_t *cr)
{
	for (const ball *b = balls; b != balls + n_balls; ++b)
		b->draw(cr);
}

void border_draw(cairo_t *cr)
{
	cairo_save(cr);
	cairo_new_path(cr);
	cairo_move_to(cr, width / 2, height / 2);
	cairo_line_to(cr, width / 2, height / 2 - border_velocity);
	cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
	cairo_set_line_width(cr, 1.0);
	cairo_stroke(cr);
	cairo_arc(cr, width / 2, height / 2 - border_velocity, 3, 0, 2 * M_PI);
	cairo_fill(cr);
	cairo_restore(cr);
}

static void balls_destroy_faces()
{
	for (ball_face *f : faces)
		if (f)
			delete (f);
	faces.clear();
}

void balls_destroy()
{
	balls_destroy_faces();
	delete[] (balls);
}

void balls_init()
{
	balls = new ball[n_balls];
	assert(balls);
	if (!track)
		balls_init_state();
	else
		balls_init_state_track();
	balls_init_faces();
}
