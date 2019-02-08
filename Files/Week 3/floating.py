import numpy as np
import matplotlib.pyplot as plt

value, float, double = abs(np.genfromtxt('floating.dat').T)

fig, ax = plt.subplots(figsize = (8, 4.5))
ax.loglog(value, float, '-', label = '32 bit')
ax.loglog(value, double, '-', label = '64 bit')

# ax.loglog(value, 100 * float / value)
# ax.loglog(value, 100 * double / value)

ax.set_xlabel('Floating point value')
ax.set_ylabel('Floating point precision')
ax.legend()

fig.savefig('floating.png', dpi = 300)

# plt.show()
