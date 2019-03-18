import sys
sys.path += ["../../.."]

import smbh
import matplotlib.pyplot as plt


results = smbh.Results("return_0_1.dat")
r = results.distance
t = results.times

plt.plot(t, r)
plt.show()
