#%%
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
#%%
N = 1000
M = 10000
data = str(subprocess.check_output(["/Users/samedi/Documents/прога/study_modelling/bins/random_walk_hist", str(N), str(M)]))[2:-2]
ser = pd.Series(map(float, data.split(" ")))
ax = ser.hist(bins=100)
plt.show()
# %%
