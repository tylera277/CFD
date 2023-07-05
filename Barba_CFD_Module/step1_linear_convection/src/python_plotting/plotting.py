
# Program: General program with the intention of plotting
#  output from separate cpp files
# Date: 04Jul2023

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Importing data from a csv exchange file
raw_data = pd.read_csv('../output_data/04Jul2023_data.csv')
print(raw_data.head())
raw_data = raw_data.dropna()
numpy_data = raw_data.to_numpy()
print("pls:", (numpy_data.transpose))


# Number of data points to plot
nx = 40
x = np.linspace(0, 2, 40)
print(x.shape)
print(numpy_data)
numpy_data = np.transpose(numpy_data)
print(numpy_data)
plt.plot(x, numpy_data[:])
plt.show()
