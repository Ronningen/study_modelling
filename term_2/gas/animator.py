import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import json

config = json.load(open("/Users/samedi/Documents/прога/study_modelling/term_2/gas/config.json"))
dots_j = config["dots"]
N = len(dots_j)
dots_type = [dots_j[i]["type"] for i in range(N)]
types = {}
for t in config["types"]:
    types[t["name"]] = t
collisions = 1000
data = (str(subprocess.check_output(["/Users/samedi/Documents/прога/study_modelling/bins/gas", str(collisions)]))[2:-1]).split('\\n')

collision = []
time = []
x = []
vx = []
energy = []

print(data[-1])

for row in data[:-1]:
    if (row.__contains__("no more collisions")):
        print("stopped")
        break
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

length = 10
fig = plt.figure(figsize=(length,7))

plt.subplot(231)
plt.plot(time, collision)
plt.grid()
plt.xlabel("time")
plt.title("collisions")

plt.subplot(232)
for r in np.array(x).T:
    plt.plot(time, r)
plt.xlabel("time")
plt.title("world lines of dots")

plt.subplot(233)
for r in np.array(energy).T:
    plt.plot(time, r)
plt.plot(time, np.array(energy).T[-1], label="full energy")
plt.grid()
plt.legend()
plt.title("energies")

ax = plt.subplot(313)
ax.set_ylim(-1,1)
ax.set_yticks([])
ln = ax.scatter([], [])
k = 1100*length/(config["box"]["r"]-config["box"]["l"])
def init():
    ax.set_xlim(config["box"]["l"], config["box"]["r"])
    ax.set_ylim(-1,1)
    ax.set_yticks([])
    ln.set_offsets(np.array([x[0], [0]*len(x[0])]).T)
    ln.set_color([types[d]["color"] for d in dots_type])
    ln.set_sizes([types[d]["radius"]*k for d in dots_type])
    return ln,

def update(t, ln, i):
    while i < len(time) - 1 and time[i + 1] <= t:
        i += 1
    data = np.array(x[i]) + np.array(vx[i]) * (t-time[i])
    ln.set_offsets(np.array([data, [0]*len(data)]).T)
    return ln,

t = time[len(time)-1]
dt = 0.1
i = 0

ani = FuncAnimation(fig, func=update, frames=np.arange(0,t,dt), fargs=(ln, i),
    init_func=init, blit=True, interval = 33)

plt.show()

    