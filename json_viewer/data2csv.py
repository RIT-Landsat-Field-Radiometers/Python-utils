import json
import csv
import sys
from os import walk
from os import path

import datetime as datetime
import matplotlib.pyplot as plt
import numpy as np

from _datetime import datetime

if __name__ == "__main__":

    # # asking for datetime range for data #
    # STRrangeStart = input("please give mm/dd/yyyy for start date")
    # print("rangeStart : " + STRrangeStart)
    # # validate start date
    # try:
    #     date = datetime.strptime(STRrangeStart, '%m/%d/%Y')
    # except ValueError:
    #     print("Invalid start date")
    # STRrangeEnd = input("please give mm/dd/yyyy for end date")
    # print("rangeEnd : " + STRrangeEnd)
    #
    # # validate end date
    # try:
    #     date = datetime.strptime(STRrangeEnd, '%m/%d/%Y')
    # except ValueError:
    #     print("Invalid end date")
    #
    # monStart, dayStart, yearStart = STRrangeStart.split('/')
    # monEnd, dayEnd, yearEnd = STRrangeEnd.split('/')

    # begin looking for files #
    mypath = "./Data/"
    # filenames = next(walk(mypath), (None, None, []))[2]  # [] if no file. this line does not work

    # START DEBUGGING WALK #
    # print(filenames)
    # filenames = list(walk(mypath))  # type: List[Tuple[str, List[str], List[str]]]    # grabs all files in entire directory
    # print(filenames)
    # print("\n")
    filenames = []
    for root, dirs, files in walk(mypath):
        for name in files:
            print(name)
            print(path.join(root, name))
            # filenames.append(name)
            filenames.append(path.join(root, name))
    print("filenames\n")
    print(filenames)
    if filenames[0] == './Data/output.csv':
        filenames = filenames[1::2]
    else:
        filenames = filenames[0::2]
    print("no rpb")
    print(filenames)
    # END DEBUGGING WALK #

    # filenames.sort(key=lambda x: int(x[1:-5]))

    first = True

    sensors = [[[] for _ in range(8)] for _ in range(2)]

    chans = [[] for _ in range(8)]
    internalTemps = [[[] for _ in range(2)] for _ in range(2)]
    hum = []
    pres = []
    rain = []
    rainTemp = []   # array to hold current json rain vals for dataTime loops
    windspd = []
    winddir = []
    airTemp = []
    events = []
    dataTime_extracted = []  # array to hold the dataStart val from json file
    dataTime_offset = []    # array to place the JSON val and add offset to
    dataTime_hr = []    # human readable datetime
    # duration = []   # array to hold duration val from json file

    for fidx, fname in enumerate(filenames):
        # with open(mypath + fname, 'r') as dfp: #original line
        with open(fname, 'r') as dfp:
            data = json.loads(dfp.read())
            print(fname)
            for sidx, sen in enumerate(data['sensors']):
                if not sen:
                    continue
                #   extract channel values from sensor data
                ichans = sen['channels']
                ichans = [x['values'] for x in ichans]
                for cidx, i in enumerate(ichans):
                    sensors[sidx][cidx].extend(i)
                # extract values from internal temp data
                iTemps = sen['internalTemps']
                iTemps = [x['values'] for x in iTemps]
                for itidx, it in enumerate(iTemps):
                    internalTemps[sidx][itidx].extend(it)
            tx = np.array(data['bmeBoard']['humidity']['values'])
            tx = tx.astype(float)
            hum.extend(tx)
            # hum.extend(data['bmeBoard']['humidity']['values'])
            pres.extend(data['bmeBoard']['pressure']['values'])
            rain.extend(data['bmeBoard']['rain'])
            windspd.extend(data['bmeBoard']['windSpeed']['values'])
            airTemp.extend(data['bmeBoard']['airTemperature']['values'])
            rainTemp.extend(data['bmeBoard']['rain'])   # gather rain vals from current json file

            tx = np.array(data['bmeBoard']['windSpeed']['values'])
            tx = tx.astype(float)
            tx = [x * (x < 50) for x in tx]
            # windspd.extend(data['bmeBoard']['windSpeed']['values'])
            winddir.extend(data['bmeBoard']['windDirection']['values'])
            # duration.extend(data['duration'])
            dataTime_extracted.extend([str(data['dataStart']['unixTime'])])  # grab initial startTime
            # print(dataTime)
            for ridx, rainStatus in enumerate(rainTemp):
                dataTime_offset.extend([ str(int(dataTime_extracted[fidx]) + (ridx+1)*30) ])    # add 30 sec per entry
                # print(dataTime)
            rainTemp = []   # reset rainTemp vals for next json file

    for tidx, timeValue in enumerate(dataTime_offset):
        dataTime_hr.extend([datetime.fromtimestamp(float(dataTime_offset[tidx])).strftime('%Y-%m-%d %H:%M:%S')])




    # writes data into single column test
    # sensorTest = [[[1, 2, 3],
    #                [4, 5, 6],
    #                [7, 8, 9],
    #                [10, 11, 12]],
    #               [[13, 14, 15],
    #                [16, 17, 18],
    #                [19, 20, 21],
    #                [22, 23, 24]]]
    # humTest = [68, 69, 70]
    # presTest = [100, 101, 102]
    # BMEdata = [hum] + [pres] + [rain] + [winddir]
    # data = sensors[0] + sensors[1] + BMEdata[0] + BMEdata[1] + BMEdata[2] + BMEdata[3]

    # getting the time column #

    BMEdata = [[hum], [pres], [rain], [windspd], [winddir], [airTemp]]
    data = [dataTime_hr] + sensors[0] + internalTemps[0] + sensors[1] + internalTemps[1] + \
        BMEdata[0] + BMEdata[1] + BMEdata[2] + BMEdata[3] + BMEdata[4] + BMEdata[5]
    headerList = ['Time'] + ['T1CA'] + ['T1CB'] + ['T1CC'] + ['T1CD'] + \
                 ['T2CA'] + ['T2CB'] + ['T2CC'] + ['T2CD'] + ['IT1'] + ['IT2'] + \
                 ['T1CA'] + ['T1CB'] + ['T1CC'] + ['T1CD'] + \
                 ['T2CA'] + ['T2CB'] + ['T2CC'] + ['T2CD'] + ['IT1'] + ['IT2'] + \
                 ['Humidity'] + ['Pressure'] + ['Rain'] + ['Wind Speed'] + ['Wind Direction'] + ['Air Temperature']
    with open("./output.csv", "w", encoding='utf8', newline='') as csvfile:
        hrd = csv.writer(csvfile)
        # hrd.writerow(["C%d" % d for d in range(1, 1 + len(data))])
        hrd.writerow(headerList)
        for r in zip(*data):
            hrd.writerow(r)

    # with open("./Data/converted.csv", "w") as ofp:
    #     ofp.write("Time,C1,C2,C3,C4,C5,C6,C7,C8\n")
    #     for toff in range(len(sensors[0])):
    #         print(toff)
    #         ofp.write(datetime.fromtimestamp(startTime + toff).strftime('%Y-%m-%d %H:%M:%S') + ",")
    #         for cidx in range(len(sensors)):
    #             ofp.write(str(sensors[cidx][toff]) + ",")
    #         ofp.write("\n")
