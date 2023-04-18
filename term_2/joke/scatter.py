import numpy as np
import matplotlib.pyplot as plt
import subprocess

plt.figure(figsize=(10,7))
plt.ylim(0,1)

x0 = 0.5
data = (str(subprocess.check_output(["/Users/samedi/Documents/прога/study_modelling/bins/joke", str(x0)]))[2:-1]).split('\\n')
for r in data[:-1]:
    r = r.split(": ")
    y = list(map(float, r[1][:-1].split(" ")))[150:]
    plt.scatter([float(r[0])]*len(y), y, s=0.03, c="black", alpha=0.2)
plt.xlabel(r"$alpha$")
plt.ylabel(r"$x_n\ 0<n<100$")
plt.title(r"$x_0=$"+str(x0))

plt.show()
