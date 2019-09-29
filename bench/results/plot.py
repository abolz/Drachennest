import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set()

names = [
    "charconv",
    "charconv-general",
    "grisu2",
    "grisu3",
    "ryu",
]

for name in names:
    data = pd.read_csv(name + '.csv')
    data = data.pivot('len(N)', 'E', 'ns')
    data = data.transpose()
    ax = sns.heatmap(data, vmin=0.0, vmax=150.0, cmap="inferno")
    plt.savefig(name + '.png')
    plt.close()
