import json
from os import walk

import matplotlib.pyplot as plt
import numpy as np

from datetime import datetime

if __name__ == "__main__":

    mypath = "./Month/"

    dirnames = next(walk(mypath), (None, None, []))[1]
    dirnames.sort(key=lambda x: int(x))

    filenames = []

    for dirName in dirnames:
        tfiles = next(walk(mypath + dirName), (None, None, []))[2]  # [] if no file
        tfiles.sort(key=lambda x: int(x[1:-5]))
        tfiles = [mypath + dirName + "/" + x for x in tfiles]
        filenames.extend(tfiles)

    first = True

    sensors = [[[] for _ in range(8)] for _ in range(2)]

    chans = [[] for _ in range(8)]
    hum = []
    pres = []
    rain = []
    windspd = []
    winddir = []
    events = []

    for fname in filenames:
        with open(fname, 'r') as dfp:
            data = json.loads(dfp.read())
            print(fname)
            for sidx, sen in enumerate(data['sensors']):
                if not sen:
                    continue
                ichans = sen['channels']
                ichans = [x['values'] for x in ichans]
                for cidx, i in enumerate(ichans):
                    sensors[sidx][cidx].extend(i)
                    # chans[cidx].extend(i)
            tx = np.array(data['bmeBoard']['humidity']['values'])
            tx = tx.astype(np.float)
            hum.extend(tx)
            # hum.extend(data['bmeBoard']['humidity']['values'])
            pres.extend(data['bmeBoard']['pressure']['values'])
            rain.extend(data['bmeBoard']['rain'])

            tx = np.array(data['bmeBoard']['windSpeed']['values'])
            tx = tx.astype(np.float)
            tx = [x * (x < 50) for x in tx]
            windspd.extend(tx)
            # windspd.extend(data['bmeBoard']['windSpeed']['values'])
            winddir.extend(data['bmeBoard']['windDirection']['values'])

    # with open("./Data/converted.csv", "w") as ofp:
    #     ofp.write("Time,C1,C2,C3,C4,C5,C6,C7,C8\n")
    #     for toff in range(len(chans[0])):
    #         ofp.write(datetime.fromtimestamp(startTime + toff).strftime('%Y-%m-%d %H:%M:%S') + ",")
    #         for cidx in range(len(chans)):
    #             ofp.write(str(chans[cidx][toff]) + ",")
    #         ofp.write("\n")

    plt.figure(1)
    plt.plot(hum)
    plt.title("Relative Humidity")
    plt.xlabel("Time(s)")
    plt.ylabel("RHUM %")

    plt.figure(2)
    plt.plot(pres)
    plt.title("Atmospheric Pressure")
    plt.xlabel("Time(s)")
    plt.ylabel("Pressure (Pa)")

    plt.figure(3)
    plt.plot(rain)
    plt.title("Rain Detected")
    plt.xlabel("Time(s)")
    plt.ylabel("Yes/No")

    plt.figure(5)
    wfig, (wax1, wax2) = plt.subplots(2, 1)
    wfig.suptitle("Wind Speed/Direction")
    wax1.plot(windspd)
    wax1.set_ylabel("Speed (m/s)")
    wax2.plot(winddir)
    wax2.set_ylabel("Degrees from North")
    wax2.set_xlabel("Time (s)")

    for sidx, sen in enumerate(sensors):
        plt.figure(6 + sidx)
        fig, ax = plt.subplots()
        for idx, c in enumerate(sen):
            ax.plot(c)
        ax.legend(['chan' + str(x) for x in range(len(sen))])
        plt.ylabel("Volts")
        plt.title("Sensor" + str(sidx))
        plt.xlabel("Seconds")
    plt.show()
    pass
