
# Program: General program with the intention of plotting
#  output from separate cpp files
# Date: 06Jul2023

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mpl_toolkits.mplot3d import Axes3D


# Importing data from a csv exchange file
raw_data = pd.read_csv('../output_data/09Jul2023_data.csv')

print("PLEASE: ", raw_data.head)

raw_data = raw_data.dropna()
numpy_data = raw_data.to_numpy()

print(numpy_data)

# Number of data points to plot
nx = 40
ny = 40

# X-points to attach the c++ program output to for plotting
# purposes
x = np.linspace(0, 2, nx)
y = np.linspace(0, 2, ny)


print("shape1: ", numpy_data.shape)
numpy_data = np.transpose(numpy_data)
print("SHAPE:", numpy_data.shape)


fig = plt.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x, y)
surf = ax.plot_surface(X, Y, numpy_data[:], cmap=plt.cm.viridis)
plt.show()

