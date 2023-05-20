CXX = g++
MAGICK = `Magick++-config --cppflags --cxxflags --ldflags --libs`
DFLAGS = -DVISUALIZE
CFLAGS = -O3 -pthread -std=c++20 -Wall -g $(MAGICK) $(DFLAGS)

SOURCES = algorithm.cpp
HOST = pologne
OBJECTS = main.o

run: main
	./main && wslview image.gif

main: $(OBJECTS) Makefile
	$(CXX) $(CFLAGS) -o main $(OBJECTS) 

main.o: algorithm.cpp Makefile
	$(CXX) -c $(CFLAGS) -o main.o algorithm.cpp

remote: $(SOURCES) Makefile
	scp -r $(SOURCES) Makefile $(HOST):~/tmp/
	ssh -t $(HOST) 'conda activate py311 && cd ~/tmp && ls && make main CXX=g++121 DFLAGS= MAGICK= && ./main'

clean:
	rm -f *.o
	rm -f main
