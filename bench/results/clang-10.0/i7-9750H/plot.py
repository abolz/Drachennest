import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

names = [
    "grisu3",
    "ryu",
    "schubfach",
    "dragonbox",
]

for name in names:
    data = pd.read_csv(name + '.csv')
    data = data.pivot('len(N)', 'E', 'ns')
    data = data.transpose()
    ax = sns.heatmap(data, vmin=0.0, vmax=40.0, cmap="inferno")
    plt.savefig(name + '.png')
    plt.close()
