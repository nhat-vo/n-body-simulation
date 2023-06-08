CXX = g++
DFLAGS = -DDEBUG
CFLAGS = -O3 -pthread -std=c++20 -Wall -g $(MAGICK) $(DFLAGS)

HOST = pologne
SOURCES = $(shell find $(shell pwd) -name '*.cpp' -o -name '*.hpp' -o -name '*.py' -o -name '.sh') Makefile benchmark.sh

OBJECTS = main.o visualizer.o barnes-hut.o single-thread.o multi-thread-1.o multi-thread-2.o barnes-hut-multi.o writer.o

main: $(OBJECTS)
	$(CXX) $(CFLAGS) -o main $(OBJECTS)

%.o: %.cpp  Makefile headers/common.hpp headers/vect.hpp headers/algorithms.hpp headers/visualizer.hpp headers/writer.hpp
	$(CXX) -c $(CFLAGS) -o $@ $<

visualize:
	make clean
	make DFLAGS='-DVISUALIZE -DDEBUG' MAGICK='`Magick++-config --cppflags --cxxflags --ldflags --libs`'

benchmark:
	make clean
	make DFLAGS=''

plot:
	make clean
	make DFLAGS='-DWRITE'


barnes-hut.o:  headers/barnes-hut.hpp


# This runs the code on the school's machine. You need to specify the host, and have set up a ssh authentication with it.
remote: $(SOURCES)
	scp -r $(SOURCES) $(HOST):~/tmp/
	ssh -t $(HOST) 'cd ~/tmp && ls && make CXX=g++121'

clean:
	rm -f *.o
	rm -f main
