#include <vector>
#include "balls.h"

extern int clist;
extern int max_rows;
extern int max_cols;

extern void cell_list_create();
extern void cell_list_fill();
extern void cell_list_destroy();
extern void cell_list_empty();
extern void cell_list_regenerate();
extern std::vector<ball*> neighbourhood(ball * p);