.PHONY: default
default: all

N:=100
S:=1
T:=0
P:=5
D:=0
C:=0

GTK_PACKAGES=gdk-pixbuf-2.0 gtk+-3.0
GTK_CFLAGS=$(shell pkg-config --cflags $(GTK_PACKAGES))
GTK_LIBS=$(shell pkg-config --libs $(GTK_PACKAGES))
 
# PROFILING_CFLAGS=-pg
CXXFLAGS=-Wall -g -std=c++11 -O2 $(PROFILING_CFLAGS) $(GTK_CFLAGS)

LIBS=$(GTK_LIBS) -lm

PROGS=balls
OBJS=balls.o c_index.o game.o gravity.o spaceship.o main.o polygon.o cell_list.o

# dependencies (gcc -MM *.cc)
balls.o: balls.cc game.h balls.h vec2d.h gravity.h polygon.h cell_list.h
polygon.o: polygon.cc balls.h vec2d.h
cell_list.o: balls.h cell_list.cc game.h
c_index.o: c_index.cc balls.h vec2d.h game.h c_index.h
game.o: game.cc game.h
gravity.o: gravity.cc gravity.h balls.h vec2d.h game.h
main.o: main.cc game.h balls.h vec2d.h c_index.h gravity.h spaceship.h polygon.h cell_list.h
spaceship.o: spaceship.cc balls.h vec2d.h game.h
stats.o: stats.cc

.PHONY: run
run: balls
	./balls fluid=1 spaceship=$S n=$N track=$T polygon=$P debugp=$D clist=$C

.PHONY: all
all: $(PROGS)

balls: $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) $(LIBS) -o $@

.PHONY: clean
clean:
	rm -f *.o $(PROGS) $(OBJS)
