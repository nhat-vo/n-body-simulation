CXX = g++
MAGICK = `Magick++-config --cppflags --cxxflags --ldflags --libs`
DFLAGS = -DVISUALIZE
CFLAGS = -O3 -pthread -std=c++20 -Wall -g $(MAGICK) $(DFLAGS)

SOURCES = algorithm.cpp vect.hpp visualizer.hpp common.hpp barnes-hut.hpp
HOST = pologne
OBJECTS = main.o visualizer.o

main: $(OBJECTS) Makefile
	$(CXX) $(CFLAGS) -o main $(OBJECTS)

run: main
	./main 1 10 && wslview image.gif

main.o: algorithm.cpp Makefile common.hpp barnes-hut.hpp
	$(CXX) -c $(CFLAGS) -o main.o algorithm.cpp

visualizer.o: visualizer.cpp Makefile common.hpp
	$(CXX) -c $(CFLAGS) -o visualizer.o visualizer.cpp


# This runs the code on the school's machine. You need to specify the host, and have set up a ssh authentication with it.
remote: $(SOURCES) Makefile
	scp -r $(SOURCES) Makefile $(HOST):~/tmp/
	ssh -t $(HOST) 'cd ~/tmp && ls && make main CXX=g++121 DFLAGS= MAGICK= && ./main'

clean:
	rm -f *.o
	rm -f main
