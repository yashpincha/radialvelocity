import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# add your file here
data = pd.read_csv()

x = np.array(data.t)
y = np.array(data.vel)
yerr = np.array(data.errvel)

# normalization reference 
x_ref = 0.5 * (x.min() + x.max())

t = np.linspace(x.min() - 5, x.max() + 5, 1000)

plt.figure(figsize=(12,6))
plt.errorbar(x, y, yerr=yerr, fmt=".k")
plt.xlabel("Time [days]")
_ = plt.ylabel("Radial velocity [m/s]")
