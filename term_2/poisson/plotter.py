#%%
import pandas as pd
import seaborn as sns

#%%
df = pd.read_csv('/Users/samedi/Documents/прога/study modelling/data.csv', header=None)
sns.heatmap(df)
# %%
df.iloc[0].plot()
# %%
