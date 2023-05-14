CXX = g++
CFLAGS = -pthread -std=c++17 -Wall -g `Magick++-config --cppflags --cxxflags --ldflags --libs`


OBJECTS = main.o

run: main
	./main && wslview image.gif

main: $(OBJECTS)
	$(CXX) $(CFLAGS) -o main $(OBJECTS) 

main.o: algorithm.cpp
	$(CXX) -c $(CFLAGS) -o main.o algorithm.cpp

clean:
	rm -f *.o
	rm -f main
