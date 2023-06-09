# N-Body Simulation

This repository contains the source code for a parallel N-Body simulator. A more descriptive description of the algorithms in this project is included in the report.
Here, we will give a brief introduction to the repository's structure.

## Instruction:
### Quickstart:
Run
```
make main DFLAGS=-DVISUALIZE
./main 3 1 1000
```
and a simulation will be saved in `image.gif`.

### Building:
The [Makefile](https://github.com/nhat-vo/n-body-simulation/blob/main/Makefile) contains several targets to build the project:
- `clean`: clean all the files generated by the build
- `main`: build the code into a `./main` executable
- `remote`: a useful build target to copy the code to a remote machine through `scp`, then compile the code. You can specify the remote machine in the `HOST` flag.
- `visualize`: build the code with visualization
- `benchmark`: build the code for benchmarking. This changes the output to be more suitable for reading from a script for plotting.

### Compile flags:
- In order to build the code with visualization, please add `-DVISUALIZE` in compilation. You can do this, for example, by uncommenting the `# DFLAGS = -DVISUALIZE` inside the Makefile.
Please note that this requires `imagemagick` installed, and the `Magick++config` script in `PATH`.
In addition, the visualization could be much slower compared to the actual computation, hence it is recommended to compile without this in benchmarking.
- The flag `-DWRITE` could be turned on to write the positions of the bodies into a ``.txt`` file. This is used for plotting the graphs for positions of the bodies and analyze the precision between different algorithms namely the leapfrog ($0$) and the Barnes-Hut ($3$) as the other three algorithms are only parallel versions of these two.
- The flag `-DDEBUG` enable verbose and more human-readable output. It is turned on by default and can be turned off by commenting the `# DFLAGS = -DDEBUG` inside the Makefile or using `make bencmark`.
Currently, the only debug information being printed when this is turned on is when a body exits the universe scope in Barnes-Hut algorithm.

### Binary usage:
When run without any argument, `./main` will print out the usage like below:
```
❯ ./main
Usage: ./main <algorithm id (default = 0)> <number of threads (default = 1)> <number of bodies (default = 100)>
Algorithms:
0: single-thread
1: multi-thread-1
2: multi-thread-2
3: barnes-hut
4: barnes-hut-multi-thread
```
You can specify the algorithm, number of threads, and number of bodies for the simulation. Please note that for `single-thread` and `barnes-hut`, the number of threads will always be 1.

### Benchmarking
The script `benchmark.sh` is used to benchmark the code with different algorithms and number of threads. One can run `bash ./benchmark.sh` to run the benchmarking script and modify the loop values in the script appropriately to benchmark the code with different number of threads and bodies.

## Code structure:
The repository is structured as follows:
- `vect.hpp` defines a 2D vector class and common operations for use throughout the code.
- `common.hpp` defines several configuration such as simulation period, step size, etc.
- `visualizer.hpp` and `visualizer.cpp` defines the code for the visualization functionalities. A brief description of how to use this is included in `visualizer.hpp`.
- `algorithms.hpp` declares the available algorithms.
In order to make the code more readable and avoid name collision, we only have the declaration here, and the implementation of each is split into a separate file appropriately named.
- `main.cpp` defines the argument parsing logic, constructs a scenario, then starts and times the simulation.
- `writer.hpp` and `writer.cpp` defines the code for writing the positions of the bodies from a simulation into a ``.txt`` file.
- `draw_graph.py` is a Python script to draw the graphs available in the report from the data in the ``.txt`` file.

Currently, the default scenario is 1 very heavy object (the sun) at the middle of the canvas, and `number_of_bodies - 1` objects (planets) revolving around it, with
their positions randomly distributed inside a square centered at the sun.
The planets' velocities are randomly generated so that (1) they moves in a counter-clockwise direction and (2) they orbits around the sun and neither escape nor crash into it.

A very simple scenario where a planet do a simple elliptical orbit around the sun is in the commented code in `main.cpp`.
In case you want to add/modify the scenario, please feel free to change the code in the `main.cpp` file.

## Gallery:
![1000 bodies spiraling around the sun](https://github.com/nhat-vo/n-body-simulation/blob/8c6317ca7d945acc58f24139f4c1a673f597bc1a/images/1000-bodies.gif)
![Bodies orbiting around the sun](https://github.com/nhat-vo/n-body-simulation/blob/8c6317ca7d945acc58f24139f4c1a673f597bc1a/image.gif)
![Binary stars](https://github.com/nhat-vo/n-body-simulation/blob/8c6317ca7d945acc58f24139f4c1a673f597bc1a/images/binary-stars-close.gif)
![Circular orbit](https://github.com/nhat-vo/n-body-simulation/blob/8c6317ca7d945acc58f24139f4c1a673f597bc1a/images/circle-orbit.gif)
