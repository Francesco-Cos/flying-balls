#include <cmath>
#include <cassert>
#include <iostream>
#include "balls.h"
#include "game.h"

int clist = 0;
int max_rows = height/rc;
int max_cols = width/rc;

std::vector<ball *> ** cell_list;

static bool is_in_bounds(int i, int j) {
    return !(i < 0 || i > max_rows - 1 || j < 0 || j > max_cols - 1);
}

void cell_list_create() {
    cell_list = new std::vector<ball *> * [max_rows];
    for (int i = 0; i < max_rows; i++) {
        cell_list[i] = new std::vector<ball *> [max_cols];
        cell_list[i]->reserve(n_balls);
    }
}

void cell_list_fill() {
    for (unsigned int i = 0; i < n_balls; ++i) {
        int row = floor(balls[i].position.x / rc);
        int column = floor(balls[i].position.y / rc);
        if (is_in_bounds(row, column))
            cell_list[row][column].push_back(&balls[i]);
    }
};

void cell_list_empty() {
    for (int i = 0; i < max_rows; i++) {
        for (int j = 0; j < max_cols; j++) {
            cell_list[i][j].clear();
        }
    }
}

void cell_list_destroy() {
    for (int i = 0; i < max_rows; i++) {
        delete[] cell_list[i];
    }
    delete[]cell_list;
}

void cell_list_regenerate() {
    cell_list_empty();
    cell_list_fill();
}



std::vector<ball*> neighbourhood(ball * p) {
    std::vector<ball*> neighbours;
    int row = floor(p->position.x / rc);
    int column = floor(p->position.y / rc);
    for (int i = -1; i < 2; i++) {
        for (int j = -1; j < 2; j++) {
            if (is_in_bounds(row + i,column + j)) {
                for (ball * ball : cell_list[row + i][column + j]) {
                    if (ball != p)
                        neighbours.push_back(ball);
                }
            }
        }
    }
    return neighbours;
};