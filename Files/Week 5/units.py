from scipy.constants import G as G_0 # m3s-2kg-1

kpc = 3.0856776e19 # in m
gyr = 60 * 60 * 24 * 365.25 * 1e9 # in s
sm = 1.98847e30 # in kg

distance = kpc ** -3
time = gyr ** 2
mass = sm # 1e5 solar masses base unit

G = G_0 * distance * time * mass
print(G, G_0)
