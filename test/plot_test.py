import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(10, 8))
tips = pd.read_csv("tips.csv")
b = sns.violinplot(x="day", y="total_bill", data=tips)
b.set_xlabel("X Label", fontsize=20)
b.set_ylabel("Y Label", fontsize=20)
b.tick_params(labelsize=15)
plt.show()
