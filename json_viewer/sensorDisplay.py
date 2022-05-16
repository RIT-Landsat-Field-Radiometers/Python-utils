import json
from os import walk

import matplotlib.pyplot as plt
import numpy as np

from datetime import datetime

if __name__ == "__main__":

    temp1 = []
    temp2 = []

    channels1 = [[], [], [], []]
    channels2 = [[], [], [], []]

    with open('./Sensordata/sensor_log.csv', 'r') as dfp:

        header = True

        while True:
            line = dfp.readline()
            if not line:
                break

            if header:
                header = False
                continue

            splitted = line.strip().split(',')

            temp1.append(float(splitted[5]))
            temp2.append(float(splitted[10]))

            # for i in range(4):
            #     channels1[i].append(splitted[i + 1])
            #     channels2[i].append(splitted[i + 6])

    plt.figure(1)
    fig, ax = plt.subplots()
    ax.plot(temp1)
    ax.plot(temp2)
    plt.show()

    pass
