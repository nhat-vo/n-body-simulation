CXX = g++
# MAGICK = `Magick++-config --cppflags --cxxflags --ldflags --libs`
MAGICK =
# DFLAGS = -DVISUALIZE
# DFLAGS = -DWRITE
# DFLAGS =
CFLAGS = -O3 -pthread -std=c++20 -Wall -g $(MAGICK) $(DFLAGS)

HOST = pologne
SOURCES = $(shell find $(SOURCE_DIR) -name '*.cpp' -o -name '*.hpp') Makefile

OBJECTS = main.o visualizer.o barnes-hut.o single-thread.o multi-thread-1.o multi-thread-2.o barnes-hut-multi.o writer.o

main: $(OBJECTS)
	$(CXX) $(CFLAGS) -o main $(OBJECTS)

%.o: %.cpp  Makefile common.hpp vect.hpp algorithms.hpp visualizer.hpp writer.hpp
	$(CXX) -c $(CFLAGS) -o $@ $<

barnes-hut.o:  barnes-hut.hpp


# This runs the code on the school's machine. You need to specify the host, and have set up a ssh authentication with it.
remote: $(SOURCES)
	scp -r $(SOURCES) $(HOST):~/tmp/
	ssh -t $(HOST) 'cd ~/tmp && ls && make main CXX=g++121 DFLAGS= MAGICK='

clean:
	rm -f *.o
	rm -f main
