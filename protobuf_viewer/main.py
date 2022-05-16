
import hourly_pb2 as hourly
from os import walk

import matplotlib.pyplot as plt


def calc_temp(Vout, Td):

    epsilon = 0.995
    F = [0.02327274, 0.99999971, 0.02953881, 0.16478234, 0.01299859, 0.01875265, 0.01527482, 0.01265653]
    F1 = [-0.02894337, -0.031758, -0.02744425, -0.00098452, -0.03816642, -0.05631432, -0.05692207, -0.07251704]
    n = [2., 1.31652025, 1.99999998, 2., 1.99999998, 2., 1.99999999, 1.99999999]

    outputs = []
    for cidx in range(len(Vout)):
        vin = (Vout[cidx] - 1.25) * 1000
        tint = 273.15
        if cidx < 4:
            tint += Td[0]
        else:
            tint += Td[1]

        outputs.append(abs(1 / epsilon * (vin / F[cidx] - (F1[cidx] * tint ** n[cidx]) + tint ** n[cidx])) ** (1 / n[cidx]))
        outputs[-1] -= 273.15

    return outputs



if __name__ == '__main__':

    mypath = "./Data/"
    filenames = next(walk(mypath), (None, None, []))[2]  # [] if no file
    filenames.sort(key=lambda x: int(x[1:-4]))

    chans = [[] for _ in range(8)]
    itemps = [[] for _ in range(2)]
    hum = []
    pres = []
    rain = []
    windspd = []
    winddir = []
    airtemp = []
    events = []

    for fname in filenames:
        protomessage = hourly.HourlyData()
        with open(mypath + fname, "rb") as fp:
            protomessage.ParseFromString(fp.read())
        for cidx in range(8):
            chans[cidx].extend(list(protomessage.sensors[0].channels[cidx].values))
        for itidx in range(2):
            itemps[itidx].extend(list(protomessage.sensors[0].internalTemps[itidx].values))

        hum.extend(list(protomessage.bmeBoard.humidity.values))
        pres.extend(list(protomessage.bmeBoard.pressure.values))
        rain.extend(list(protomessage.bmeBoard.rain))
        windspd.extend(list(protomessage.bmeBoard.windSpeed.values))
        winddir.extend(list(protomessage.bmeBoard.windDirection.values))
        airtemp.extend(list(protomessage.bmeBoard.airTemperature.values))

        # if len(hum) > 1600:
        #     break

    print(airtemp)

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
    plt.plot(airtemp)
    plt.title("Air Temperature")
    plt.xlabel("Time(s)")
    plt.ylabel("Temperature (C)")

    plt.figure(4)
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

    plt.figure(6)
    fig, ax = plt.subplots()
    leg = []
    for idx, c in enumerate(chans):
        # if idx == 7:
        #     continue
        leg.append('chan' + str(idx))
        ax.plot(c)
    ax.legend(leg)
    plt.ylabel("Volts")
    plt.title("Sensor0")
    plt.xlabel("Seconds")

    plt.figure(7)
    fig, ax = plt.subplots()
    leg = []
    for idx, c in enumerate(itemps):
        leg.append('thermopile' + str(idx))
        ax.plot(c)
    ax.legend(leg)
    plt.ylabel("Temperature (C)")
    plt.title("Sensor0 Internal Temps")
    plt.xlabel("Seconds")

    combined = list(zip(*chans))
    itempsc = list(zip(*itemps))

    zipped = list(zip(combined, itempsc))

    convtemps = []

    for pair in zipped:
        convtemps.append(calc_temp(pair[0], pair[1]))

    res = list(zip(*convtemps))

    plt.figure(8)
    fig, ax = plt.subplots()
    leg = []
    for idx, chan in enumerate(res):
        # if idx == 7:
        #     continue
        # if idx == 0:
        #     continue
        ax.plot(chan)
        leg.append('chan' + str(idx))
    ax.legend(leg)
    plt.ylabel("Degrees C")
    plt.title("Converted Values")
    plt.xlabel("Seconds")

    plt.show()

    with open("./out.csv", "w") as fp:
        fp.write("Time, Pressure, Humidity, Wind Speed, Wind Direction, Rain, Air Temperature, Internal Temp0, Internal Temp1, ")
        fp.write(", ".join(["ch" + str(i) for i in range(8)]))
        fp.write("\n")

        for idx in range(len(airtemp)):
            fp.write(", ".join([ str(val) for val in [idx*30, pres[idx], hum[idx], windspd[idx], winddir[idx],
                                                      rain[idx], airtemp[idx], itemps[0][idx], itemps[1][idx],
                                                      *[ch[idx] for ch in chans]  ]]))
            fp.write("\n")


    pass