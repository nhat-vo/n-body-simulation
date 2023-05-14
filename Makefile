CXX = g++
CFLAGS = -pthread -std=c++17 -Wall -g

OBJECTS = main.o

run: main
	./main

main: $(OBJECTS)
	$(CXX) $(CFLAGS) -o main $(OBJECTS) 

main.o: algorithm.cpp
	$(CXX) -c $(CFLAGS) -o main.o algorithm.cpp

clean:
	rm -f *.o
	rm -f main
