all: lsh cube cluster

lsh: lsh.o help_functions.o calculations.o calculations_lsh.o
	g++ lsh.o help_functions.o calculations.o calculations_lsh.o -o lsh

cube: cube.o help_functions.o calculations.o calculations_cube.o
	g++ cube.o help_functions.o calculations.o calculations_cube.o -o cube

cluster: cluster.o help_functions.o calculations.o calculations_cube.o calculations_lsh.o calculations_cluster.o
	g++ cluster.o help_functions.o calculations.o calculations_cube.o calculations_lsh.o calculations_cluster.o -o cluster

lsh.o: lsh.cpp
	g++ -c lsh.cpp

cube.o: cube.cpp
	g++ -c cube.cpp

cluster.o: cluster.cpp
	g++ -c cluster.cpp

help_functions.o: help_functions.cpp
	g++ -c help_functions.cpp

calculations.o: calculations.cpp
	g++ -c calculations.cpp
	
calculations_lsh.o: calculations_lsh.cpp
	g++ -c calculations_lsh.cpp

calculations_cube.o: calculations_cube.cpp
	g++ -c calculations_cube.cpp

calculations_cluster.o: calculations_cluster.cpp
	g++ -c calculations_cluster.cpp
	
clean:
	rm -f lsh cube cluster help_functions calculations calculations_lsh calculations_cube calculations_cluster lsh.o cube.o cluster.o help_functions.o calculations.o calculations_lsh.o calculations_cube.o calculations_cluster.o
