import json
from os import walk

import matplotlib.pyplot as plt
import numpy as np

from datetime import datetime

if __name__ == "__main__":

    target = []
    rtd = []
    internals = [[], []]
    externals = [[] for _ in range(8)]
    bbtemp = []

    with open("./Data2/sensor_log.csv", "r") as ifp:
        ifp.readline()
        while True:
            line = ifp.readline()
            if not line:
                break
            data = line.split(",")
            target.append(float(data[2]))

            rtd.append(float(data[3]))
            bbtemp.append(float(data[1]))

            internals[0].append(float(data[8]))
            internals[1].append(float(data[13]))

            externals[0].append(float(data[4]))
            externals[1].append(float(data[5]))
            externals[2].append(float(data[6]))
            externals[3].append(float(data[7]))

            externals[4].append(float(data[9]))
            externals[5].append(float(data[10]))
            externals[6].append(float(data[11]))
            externals[7].append(float(data[12]))

    print("Done")

    plt.figure(1)
    plt.plot(target)
    plt.plot(rtd)
    plt.plot(internals[0])
    plt.plot(internals[1])
    plt.title("Temps")
    plt.xlabel("Time(s)")
    plt.ylabel("Deg *C")
    plt.legend(["Target", "RTD", "INT1", "INT2"])

    plt.figure(2)
    plt.plot(bbtemp)
    plt.title("Black Body")
    plt.xlabel("Time(s)")
    plt.ylabel("Deg *C")

    plt.figure(3)
    plt.title("Sensor 1")
    plt.xlabel("Time(s)")
    plt.ylabel("Volts")
    for i in range(4):
        plt.plot(externals[i])
    plt.legend(["chan0", "chan1", "chan2", "chan3"])

    plt.figure(4)
    plt.title("Sensor 2")
    plt.xlabel("Time(s)")
    plt.ylabel("Volts")
    for i in range(4):
        plt.plot(externals[i + 4])
    plt.legend(["chan0", "chan1", "chan2", "chan3"])

    plt.show()

    pass