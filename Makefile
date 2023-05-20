CXX = g++
CFLAGS = -O3 -pthread -std=c++20 -Wall -g `Magick++-config --cppflags --cxxflags --ldflags --libs`

SOURCES = algorithm.cpp
HOST = pologne
OBJECTS = main.o

run: main
	./main && wslview image.gif

main: $(OBJECTS)
	$(CXX) $(CFLAGS) -o main $(OBJECTS) 

main.o: algorithm.cpp
	$(CXX) -c $(CFLAGS) -o main.o algorithm.cpp

remote: $(SOURCES) Makefile
	scp -r $(SOURCES) Makefile $(HOST):~/tmp/
	ssh -t $(HOST) 'conda activate py311 && cd ~/tmp && ls && make main CXX=g++121 && ./main'
	rm -f image.png
	scp $(HOST):~/tmp/image.png $(SOURCE_DIR)
	wslview image.png

clean:
	rm -f *.o
	rm -f main
