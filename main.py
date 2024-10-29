tmax=18 # days
vmax=1e4 # km/s
r = vmax * tmax * 3600 * 24

z=1200 # km/s
d_z = z/ 70 # Mpc
d_z = d_z * 3.09e13  # Mpc to km

theta = r/d_z
print(theta)

d = 176 / (theta * 1e6)
print(d)
