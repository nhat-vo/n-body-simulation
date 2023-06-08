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

    plt.figure(figsize=(6, 6))
    plt.plot(error, label="error")
    plt.xlabel("Time")
    plt.ylabel("Error")
    plt.legend()
    plt.savefig("graphs/error.png")

    # plot the graph
    plt.figure(figsize=(6, 6))
    plt.plot(distance_single, label="single")
    plt.plot(distance_barnes, label="barnes")
    plt.xlabel("Time")
    plt.ylabel("Distance")
    plt.legend()
    plt.savefig("graphs/distance.png")

    # plot the graph
    plt.figure(figsize=(6, 6))
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
        plt.figure(figsize=(6, 6))
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

def thread_algo_graph():
    data = np.loadtxt("thread_algo.txt", dtype=float, delimiter=" ")
    multi_algo_data = np.array([data[:4, 1:], data[4:8, 1:], data[8:12, 1:]])

    # print(multi_algo_data[1, ])
    # return
    threads = [1, 2, 4, 8]

    algo_names = {0: "multi-1", 1: "multi-2", 2: "barnes-multi"}
    for algo in range(3):
        plt.figure(figsize=(6, 6))
        plt.plot(threads, multi_algo_data[algo, 0, :], label="5 bodies")
        plt.plot(threads, multi_algo_data[algo, 1, :], label="50 bodies")
        plt.plot(threads, multi_algo_data[algo, 2, :], label="500 bodies")
        if algo == 2:
            plt.plot(threads, multi_algo_data[algo, 3, :], label="5000 bodies")
        plt.xlabel("Threads")
        plt.ylabel("Time")
        plt.title("{} algorithm".format(algo_names[algo]))
        # plt.yscale("log")
        plt.legend()
        plt.savefig("graphs/thread_algo_{}.png".format(algo_names[algo]))

if __name__ == "__main__":
    # thread_algo_graph()
    simulation_time_graphs()
    # error_graphs()