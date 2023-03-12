#%%
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
#%%
N = 1000
M = 20000
data = str(subprocess.check_output(["/Users/samedi/Documents/прога/study_modelling/bins/random_walk", str(N), str(M)]))[2:-2]
data = data.split("|")
steps = range(1, N+1)
offsets = list(map(float, data[0][:-1].split(',')))
offsets_abs = list(map(float, data[1][:-1].split(',')))
offsets_2 = list(map(float, data[2][:-1].split(',')))

plt.grid(True)
plt.plot(steps, offsets, label="смещение")
plt.plot(steps, offsets_abs, label="модуль смещения")
plt.plot(steps, offsets_2, label="квадрат смешения")

# plt.plot(steps, np.array(offsets_abs)/(np.sqrt(2/np.pi)*np.sqrt(steps))-1, label="отклонение модуля")
# plt.plot(steps, np.array(offsets_2)/steps-1, label="отклонение квадрата")
# plt.legend()
plt.show()
# %%
