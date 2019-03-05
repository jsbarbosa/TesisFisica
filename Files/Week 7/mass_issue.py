import numpy as np
import matplotlib.pyplot as plt

r = 2
a = 1
b = 2

theta = np.linspace(0, 2 * np.pi)

x = (r / a) * np.cos(theta)
y = (r / b) * np.sin(theta)
x2 = (r / b) * np.cos(theta)
y1 = (r / a) * np.sin(theta)

fig, ax = plt.subplots()

ax.plot(x, y, c = 'k')
ax.plot(x, y1)
ax.plot(x2, y)

ax.plot([0, 0], [1.1 * r / a, -1.1 * r / a], c = 'k', lw = 0.5)
ax.plot([1.1 * r / a, -1.1 * r / a], [0, 0], c = 'k', lw = 0.5)


ax.scatter([r / a], [0])
ax.scatter([0], [r / b])

ax.axis('off')
ax.axis('equal')

fig.tight_layout()

fig.savefig("triaxial_mass_issue.png", dpi = 300)
