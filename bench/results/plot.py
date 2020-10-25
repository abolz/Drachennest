import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# uniform (0,1)
data = [
    # gr2*     gr2       gr3        ryu        sf        db       db*      std
    [ 25.0,    29.0,     25.1,      25.1,      20.4,     25.4,    16.6,    49.3],
    [ 25.8,    29.9,     25.3,      25.2,      20.5,     27.0,    15.8,    49.3],
    [ 24.3,    28.5,     23.9,      26.7,      20.8,     25.9,    16.1,    47.1],
    [ 22.5,    27.1,     22.5,      25.3,      21.9,     26.4,    17.5,    45.4],
    [ 21.7,    26.5,     23.7,      25.9,      22.0,     26.4,    17.6,    44.7],
    [ 20.9,    25.8,     25.1,      25.2,      21.6,     25.9,    17.1,    42.4],
    [ 21.2,    26.2,     27.0,      26.0,      21.6,     26.1,    16.9,    42.4],
    [ 20.7,    25.1,     28.1,      25.4,      22.8,     27.3,    18.3,    40.8],
    [ 21.5,    27.2,     30.6,      26.4,      23.6,     28.9,    18.6,    38.4],
    [ 22.7,    27.2,     33.1,      26.6,      23.0,     30.1,    18.5,    39.4],
    [ 21.8,    26.5,     32.5,      25.5,      23.4,     29.0,    18.5,    37.6],
    [ 21.7,    26.3,     34.5,      25.2,      23.5,     29.0,    19.2,    35.7],
    [ 23.8,    28.2,     36.2,      25.4,      24.4,     27.9,    19.6,    34.2],
    [ 23.5,    27.6,     37.5,      24.5,      23.7,     28.8,    19.3,    33.1],
    [ 24.4,    28.3,     39.6,      25.3,      24.2,     28.3,    19.7,    32.6],
    [ 23.6,    28.2,     41.2,      24.7,      25.2,     27.3,    21.4,    30.4],
    [ 24.6,    29.8,     45.8,      25.5,      25.7,     27.7,    23.3,    30.7],
]

index = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
df = pd.DataFrame(data, index, ["Grisu2*", "Grisu2**", "Grisu3", "Ryu", "Schubfach", "Dragonbox", "Dragonbox*", "std::to_chars"])

ax = sns.lineplot(data=df, dashes=False)
plt.savefig('bench' + '.png')
plt.close()

# names = [
#     "charconv",
#     "charconv-general",
#     "grisu2",
#     "grisu3",
#     "ryu",
# ]

# for name in names:
#     data = pd.read_csv(name + '.csv')
#     data = data.pivot('len(N)', 'E', 'ns')
#     data = data.transpose()
#     ax = sns.heatmap(data, vmin=0.0, vmax=150.0, cmap="inferno")
#     plt.savefig(name + '.png')
#     plt.close()
