import subprocess
import numpy as np
import matplotlib.pyplot as plt
import json

data = str(subprocess.check_output(f"./bins/solver"))[2:-1]
conf = json.load(open("term_1/solver/config.json"))
head = conf["head"]
if head["do_log"]:
    el_sep = head["el_sep"]
    zone_sep = head["zone_sep"]
    row_sep = "\\n" # head["row_sep"] # TODO: automatic control of .split("\n")
    run_sep = head["run_sep"]

    exps = data.split(run_sep)[:-1]
    # TODO: add plotting of analytical solution

    separete = head["show_separate"]
    if separete:
        plotrows = len(exps)
    else:
        plotrows = 1
    plotrow = 0

    for exp in exps:
        X = []
        Y1 = []
        Y2 = []
        I = []

        for row_ in exp.split(row_sep)[1:-1]: # TODO: add first row of first experiment
            zones = list(map(lambda l: str.strip(l, el_sep), row_.split(zone_sep)))
            first = list(map(float, zones[0].split(el_sep)))
            X.append(first[0])
            Y1.append(first[1]) 
            Y2.append(first[2]) 
            if len(zones) > 1:
                I.append(float(zones[1].split(el_sep)[0]))
        
        plt.subplot(plotrows, 2, plotrow * 2 + 1)
        plt.plot(X, Y1)
        plt.grid(True)
        plt.subplot(plotrows, 2, plotrow * 2 + 2)
        plt.plot(X, I)
        # plt.plot(Y1, Y2)
        # plt.grid(True)
        if separete:
            plotrow += 1
    
plt.show()
