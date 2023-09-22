# This python file needs to loop through testing folders and creat ssingle csv file for each method

import os
import numpy as np
import matplotlib.pyplot as plt

def main():
    num_trajecs = 100
    methods = ["baseline", "SI2"]
    labels = ["optTime", "costReduction", "derivsTime", "qpTime"]

    data = np.zeros((num_trajecs, len(labels) * len(methods)))



    for i in range(len(methods)):
        for j in range(num_trajecs):

            # root = "home/davidrussell/cito_ws/src/cito/src"
            file_name = "../testingData/data/" + methods[i] + "/" + str(j) + ".csv"
            tempData = np.genfromtxt(file_name, delimiter=",")

            data[j, i * len(labels) : (i + 1) * len(labels)] = tempData

    print(data)


    optTimes = np.zeros((num_trajecs, len(methods)))
    costReductions = np.zeros((num_trajecs, len(methods)))

    for i in range(len(methods)):
        optTimes[:, i] = data[:, i * len(labels)]
        costReductions[:, i] = data[:, i * len(labels) + 1]

    fig, axes = plt.subplots(2, 1, figsize = (18,8))
    boxPlotTitle = "Average percentage calculated derivatives against interpolation methods " + "panda_pushing_clutter"
    yAxisLabel = "Average percentage calculate derivatives"
    orange = "#edb83b"
    bp3 = box_plot(avgPercentageDerivs, orange, yAxisLabel, axes[0], labels, False)

    boxPlotTitle = "average time getting derivatives against interpolation methods " + "panda_pushing_clutter"
    yAxisLabel = "Average time getting derivatives (ms)"
    orange = "#edb83b"
    bp4 = box_plot(avgTimeGettingDerivs, orange, yAxisLabel, axes[1], labels)
    fig.suptitle(taskName + " - derivative information", fontsize=16)

    # Save data to new csv file
    np.savetxt("../testingData/data.csv", data, delimiter=",")

def box_plot(data, fill_color, yAxisTitle, ax, labels, logyAxis = False, baseline_yLine = False):
    normalPosterColour = "#103755"
    highlightPosterColor = "#EEF30D"


    bp = ax.boxplot(data, patch_artist=True, meanprops={"marker":"s","markerfacecolor":highlightPosterColor, "markeredgecolor":highlightPosterColor}, showmeans=True, showfliers=False)
    if logyAxis:
        ax.set_yscale('log')
    black = "#1c1b1a"

    for element in ['medians']:
        plt.setp(bp[element], color=black)

    for element in ['means']:
        plt.setp(bp[element], color=highlightPosterColor)
        # element.set(color="#a808c4")

    # for bpMean in bp['means']:
    #     bpMean.set(color="#a808c4")

    # for element in ['boxes', 'whiskers', 'fliers', 'means', 'medians', 'caps']:
    #     plt.setp(bp[element], color=black)


    dynamicColor = "#ff8400"
    baselineColor = "#0057c9"
    setIntervalColour = "#9d00c9"
    backgroundColour = "#d6d6d6"

    
    index = 0
    for patch in bp['boxes']:
        patch.set(facecolor=normalPosterColour)

    labelSize = 11

    ax.set_ylabel(yAxisTitle, fontsize=labelSize)

    if(baseline_yLine):
        ax.axhline(y=0, color=baselineColor, linewidth=1.5, alpha=0.5)

    xticks = []
    for i in range(len(labels)):
        xticks.append(i + 1)

    ax.set_xticks(xticks)
    ax.set_xticklabels(labels, fontsize=11)
        
    return bp 
    


if __name__ == "__main__":
    main()
