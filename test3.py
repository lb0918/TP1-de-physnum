import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import moyal
from scipy.stats import skew
from scipy import integrate
import timeit

moy = moyal(loc = 150, scale = 4).rvs(size = 10000)
plt.figure(figsize = (15, 15))
plt.hist(moy, bins = 100)
plt.show()
a = skew(moy)
print('''L'asymétrie de la distribution d'énergies de moyal est de '''+str(a)+'.')