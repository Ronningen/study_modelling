import subprocess
import pandas as pd
import matplotlib.pyplot as plt

N = 10000
mean = 0
sigma = 2
data = str(subprocess.check_output(["/Users/samedi/Documents/прога/study_modelling/bins/rand_lab", str(3), str(N), str(mean), str(sigma)]))[2:-2]
# ser = pd.Series(map(float, data.split(" ")))
# ax = ser.hist(bins=100)
processed = data.split("bins: ")[1].split("weights: ")
bins = list(map(float, processed[0][:-3].split(" ")))
weights = list(map(float, processed[1].split(" ")))

plt.hist(bins[:-1], bins, weights=weights, density=True)
plt.show()