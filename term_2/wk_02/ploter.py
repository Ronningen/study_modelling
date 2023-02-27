import subprocess
import numpy as np
import matplotlib.pyplot as plt

K = 20
data = str(subprocess.check_output(["/Users/samedi/Documents/прога/study_modelling/bins/task", str(K)]))[2:-2]
data = data.split("|")
data1 = list(map(float, data[0].split()))
data2 = list(map(float, data[1].split()))

plt.grid(True)
plt.ylim(0,1)
plt.plot(range(1, len(data2)+1), data2)
plt.plot(range(1, len(data1)+1), data1)
plt.show()