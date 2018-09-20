CC      = g++
CFLAGS  = -Wall -O3 -std=c++0x -lgmp
RM      = rm -f

default: 

weight: weight.cpp
	$(CC) weight.cpp tools.cpp -o weight.o $(CFLAGS)

distance: distance.cpp
	$(CC) distance.cpp tools.cpp -o distance.o $(CFLAGS)


groebner: groebner.cpp
	$(CC) groebner.cpp tools.cpp -o groebner.o $(CFLAGS)

list: generate_list.cpp
	$(CC) generate_list.cpp -o list.o $(CFLAGS)

clean: 
	rm *.o
