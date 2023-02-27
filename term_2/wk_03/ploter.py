import subprocess
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

N = 20
data = str(subprocess.check_output(["/Users/samedi/Documents/прога/study_modelling/bins/random_walk", str(N), str(10000)]))[2:-2]
print(data)
data = data.split("|")
offsets = list(map(float, data[0][:-1].split(',')))
offsets_abs = list(map(float, data[1][:-1].split(',')))
offsets_2 = list(map(float, data[2][:-1].split(',')))

plt.grid(True)
plt.plot(range(1, N+1), offsets, label="смещение")
plt.plot(range(1, N+1), offsets_abs, label="модуль смещения")
plt.plot(range(1, N+1), offsets_2, label="квадрат смешения")
plt.legend()
plt.show()