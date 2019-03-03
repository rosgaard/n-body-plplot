all:	
	clear; \
	g++ -std=c++11 -o main main.cpp -I/usr/include/plplot -L/usr/lib64 -lplplotcxx -lplplot
