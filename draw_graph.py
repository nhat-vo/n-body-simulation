# draw grpahs from the output.txt file using matplotlib

import matplotlib.pyplot as plt
import numpy as np


def parse(f):
    # read the output.txt file
    f = open(f, "r")
    lines = f.readlines()
    f.close()

    data = []

    for line in lines:
        # split the line with space
        line = line.split()
        x, y = 0, 0
        for i in range(0, len(line), 2):
            # convert the string to float
            x += float(line[i])
            y += float(line[i+1])
        # calculate the average
        x /= len(line) // 2
        y /= len(line) // 2
        data.append([x, y])

    return np.array(data)


if __name__ == "__main__":
    # parse the output.txt file
    data_single = parse("single.txt")
    data_barnes = parse("barnes.txt")

    single_sq = np.sum((data_single)**2, axis=1)
    barnes_sq = np.sum((data_barnes)**2, axis=1)

    distance_single = np.sqrt(single_sq)
    distance_barnes = np.sqrt(barnes_sq)

    error = np.sum((data_single - data_barnes)**2, axis=1) / len(data_single)
    error_max = np.max((data_single - data_barnes)**2, axis=1)

    plt.figure(figsize=(10, 10))
    plt.plot(error, label="error")
    plt.xlabel("Time")
    plt.ylabel("Error")
    plt.legend()
    plt.savefig("error.png")

    # plot the graph
    plt.figure(figsize=(10, 10))
    plt.plot(distance_single, label="single")
    plt.plot(distance_barnes, label="barnes")
    plt.xlabel("Time")
    plt.ylabel("Distance")
    plt.legend()
    plt.savefig("distance.png")

    # plot the graph
    plt.figure(figsize=(10, 10))
    plt.plot(error_max, label="error_max")
    plt.xlabel("Time")
    plt.ylabel("Error")
    plt.legend()
    plt.savefig("error_max.png")