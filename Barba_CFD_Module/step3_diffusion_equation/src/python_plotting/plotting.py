
# Program: General program with the intention of plotting
#  output from separate cpp files
# Date: 06Jul2023

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Importing data from a csv exchange file
raw_data = pd.read_csv('../output_data/06Jul2023_data.csv')

raw_data = raw_data.dropna()
print(raw_data.head)
numpy_data = raw_data.to_numpy()


# Number of data points to plot
nx = 40

# X-points to attach the c++ program output to for plotting
# purposes
x = np.linspace(0, 2, 38)

print(numpy_data)
numpy_data = np.transpose(numpy_data)
print("SHAPE:", numpy_data.shape)

plt.plot(x, numpy_data[:])
plt.show()
