import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import json

config = json.load(open("/Users/samedi/Documents/прога/study_modelling/term_2/gas/config.json"))
dots_j = config["dots"]
N = len(dots_j)
dots_type = [dots_j[i]["type"] for i in range(N)]
collisions = 100
data = (str(subprocess.check_output(["/Users/samedi/Documents/прога/study_modelling/bins/gas", str(collisions)]))[2:-1]).split('\\n')

collision = []
time = []
x = []
vx = []
energy = []

print(data[-1])

for row in data[:-1]:
    row = row.split("collisions: ")[1]

    tmp = row.split(" time: ")
    collision.append(int(tmp[0]))
    row = tmp[1]

    tmp = row.split(" state: ")
    time.append(float(tmp[0]))
    row = tmp[1]

    tmp = row.split(" energy: ")
    state = tmp[0].split(",")[:-1]
    x.append([float(state[2*i]) for i in range(N)])
    vx.append([float(state[2*i+1]) for i in range(N)])
    row = tmp[1]

    energy.append(list(map(float, row.split(",")[:-1])))

abstime = time[0]
for i in range(1, len(time)):
    abstime += time[i]
    time[i] = abstime

plt.subplot(131)
plt.plot(time, collision)
# plt.yticks(np.arange(0, collision[-1]+1, step=1))
plt.grid()
plt.xlabel("time")
plt.title("collisions")

plt.subplot(132)
for r in np.array(x).T:
    plt.plot(time, r)
plt.xlabel("time")
plt.title("world lines of dots")

plt.subplot(133)
for r in np.array(energy).T:
    plt.plot(time, r)
plt.grid()
plt.title("energies")

plt.show()

    