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


def error_graphs():
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
    plt.savefig("graphs/error.png")

    # plot the graph
    plt.figure(figsize=(10, 10))
    plt.plot(distance_single, label="single")
    plt.plot(distance_barnes, label="barnes")
    plt.xlabel("Time")
    plt.ylabel("Distance")
    plt.legend()
    plt.savefig("graphs/distance.png")

    # plot the graph
    plt.figure(figsize=(10, 10))
    plt.plot(error_max, label="error_max")
    plt.xlabel("Time")
    plt.ylabel("Error")
    plt.legend()
    plt.savefig("graphs/error_max.png")

def simulation_time_graphs():
    data = np.loadtxt("simulation_time.txt", dtype=float, delimiter=" ")
    algo_data = np.array([data[:4, 1:], data[4:8, 1:], data[8:12, 1:], data[12:16, 1:], data[16:20, 1:]])

    # print(algo_data[1, 3, :])
    threads = [1, 2, 4, 8]

    body_pos = {0: 5, 1: 50, 2: 500, 3: 5000}
    for body in range(4):
        plt.figure(figsize=(10, 10))
        plt.plot(threads, algo_data[0, body, :], label="single")
        plt.plot(threads, algo_data[1, body, :], label="multi-1")
        plt.plot(threads, algo_data[2, body, :], label="multi-2")
        plt.plot(threads, algo_data[3, body, :], label="barnes")
        plt.plot(threads, algo_data[4, body, :], label="barnes-multi")
        plt.xlabel("Threads")
        plt.ylabel("Time")
        plt.title("Simulation Time for {} bodies".format(body_pos[body]))
        # plt.yscale("log")
        plt.legend()
        plt.savefig("graphs/simulation_time_{}.png".format(body_pos[body]))



if __name__ == "__main__":
    simulation_time_graphs()