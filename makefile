CC      = g++
CFLAGS  = -Wall -O3 -std=c++0x -lgmp
RM      = rm -f

default: all

all: guess

guess: main.cpp
	$(CC) main.cpp -o o.o $(CFLAGS)

list: generate_list.cpp
	$(CC) generate_list.cpp -o list.o $(CFLAGS)

