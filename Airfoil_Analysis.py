import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error

import matplotlib.pyplot as plt
import seaborn as sns

# Extracting data from the csv file (read csv file into pandas dataframe)
airfoil_data = pd.read_csv("AirfoilSelfNoise.csv")

# See the datatype of the input variables (frequency, angle of attack, chord length, free-stream velocity, displacement thickness)
airfoil_data.info()

# Visualize the pandas dataframe with the variables
print(airfoil_data)

# Statistical view on the data (includes percentiles, standard deviation, mean/min/max values and count/number of values)
print(airfoil_data.describe())

# Evaluate the correlation between the independent variables and dependent variable
print(airfoil_data.corr())

# Plot of each variable
sns.distplot(a=airfoil_data['f'], hist=False)
plt.show()
sns.distplot(a=airfoil_data['alpha'], hist=False)
plt.show()
sns.distplot(a=airfoil_data['c'], hist=False)
plt.show()
sns.distplot(a=airfoil_data['U_infinity'], hist=False)
plt.show()
sns.distplot(a=airfoil_data['delta'], hist=False)
plt.show()
sns.distplot(a=airfoil_data['SSPL'], hist=False)
plt.show()
