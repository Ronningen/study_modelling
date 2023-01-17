import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
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

    for exp, run in zip(exps, conf["runs"]):
        X = []
        Y1 = []
        Y2 = []
        I = []

        for row_ in exp.split(row_sep)[1:-1]: # TODO: add first row of first experiment
            zones = list(map(lambda l: str.strip(l, el_sep), row_.split(zone_sep)))
            first = list(map(float, zones[0].split(el_sep)))
            # print(first)
            X.append(first[0])
            Y1.append(first[1]) 
            Y2.append(first[2]) 
            if len(zones) > 1:
                I.append(float(zones[1].split(el_sep)[0]))

        match run["problem"]["type"]:
            case "dual_pendulum":
                plt.figure(figsize=(4,4))
                plt.xlim(-2.1, 2.1)
                plt.ylim(-2.1, 2.1)
                plt.plot(l1*np.sin(Y1), -l1*np.cos(Y1), lw=0.5)
                plt.plot(l1*np.sin(Y1)+l2*np.sin(Y2), -l1*np.cos(Y1)-l2*np.cos(Y2))
                plt.grid(True)
                l1 = run["problem"]["l1"]
                l2 = run["problem"]["l2"]
                fig, ax = plt.subplots()
                ln, = plt.plot([], [], '-o', markersize=10, color="black")
                def init():
                    ax.set_xlim(-2.1, 2.1)
                    ax.set_ylim(-2.1, 2.1)
                    fig.set_size_inches(4, 4)
                    ax.grid(True)
                    return ln,
                def update(frame):
                    ln.set_data(
                        [0, l1*np.sin(Y1[frame]), l1*np.sin(Y1[frame])+l2*np.sin(Y2[frame])], 
                        [0, -l1*np.cos(Y1[frame]), -l1*np.cos(Y1[frame])-l2*np.cos(Y2[frame])])
                    return ln,
                ani = FuncAnimation(fig, update, frames=np.arange(0,len(Y1)),
                    init_func=init, blit=True, interval = 1)
            case "ring_spring":
                T1 = (np.array(Y1)+np.array(Y2))
                T2 = (np.array(Y1)-np.array(Y2))

                fig, ax = plt.subplots()
                ln, = plt.plot([], [], '-o', markersize=10, color="black")
                def init():
                    ax.set_xlim(-1.1, 1.1)
                    ax.set_ylim(-1.1, 1.1)
                    fig.set_size_inches(4, 4)
                    ax.grid(True)
                    return ln,
                def update(frame):
                    ln.set_data(
                        [np.cos(T1[frame]), -np.cos(T2[frame])], 
                        [np.sin(T1[frame]), np.sin(T2[frame])])
                    return ln,
                ani = FuncAnimation(fig, update, frames=np.arange(0,len(Y1)),
                    init_func=init, blit=True, interval = 1)
            case "kapitza":
                w = run["problem"]["w"]
                a = run["problem"]["a"]
                l = run["problem"]["l"]

                fig, (ax0, ax) = plt.subplots(1,2)
                ln, = ax.plot([], [], '-o', markersize=10, color="black")
                def init():
                    ax.set_xlim(-1.1*l, 1.1*l)
                    ax.set_ylim(-1.1*l, 1.1*l)
                    ax0.set_xlim(-1.1*l, 1.1*l)
                    ax0.set_ylim(-1.1*l, 1.1*l)
                    fig.set_size_inches(12, 6)
                    ax0.plot(l*np.sin(Y1), a*np.cos(w*np.array(X)) - l*np.cos(np.array(Y1)), lw=0.1)
                    ax0.grid(True)
                    ax.grid(True)
                    return ln,
                def update(frame):
                    frame*=10
                    ln.set_data(
                        [0, 0, l*np.sin(Y1[frame])], 
                        [0, a*np.cos(w*X[frame]), a*np.cos(w*X[frame]) - l*np.cos(Y1[frame])])
                    return ln,
                ani = FuncAnimation(fig, update, frames=np.arange(0,len(Y1)),
                    init_func=init, blit=True, interval = 1)
            case _:
                plt.subplot(plotrows, 2, plotrow * 2 + 1)
                plt.plot(X, Y1)
                plt.grid(True)
                plt.subplot(plotrows, 2, plotrow * 2 + 2)
                plt.plot(X, I)
                plt.plot(Y1, Y2)
            
        if separete:
            plotrow += 1
    plt.show()
