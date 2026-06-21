import numpy as np
import matplotlib.pyplot as plt

mu = 398600.44e9  # Earth's gravitational parameter (m^3/s^2)

sim_time = 2* 3600*1.5   # seconds
dt = 100          # time step
n_steps = int(sim_time / dt)

# Arrays to store position and velocity
r = np.zeros((n_steps, 3))
v = np.zeros((n_steps, 3))

# Initial conditions
r[0] = [6993000, 0, 0]        # meters
v[0] = [0, 5341.2, 5341.2]    # m/s


# Orbit simulation
for i in range(n_steps - 1):

    # acceleration due to gravity
    a = -mu * r[i, :] / (np.linalg.norm(r[i, :]) ** 3)

    # update velocity
    
    v[i+1] = v[i] + a * dt

    # update position
    r[i+1] = r[i] + v[i+1] * dt


# 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot(
    r[:, 0],
    r[:, 1],
    r[:, 2]
)

# Show full numbers, not 1e7
ax.ticklabel_format(style='plain')

# Labels
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.set_zlabel('Z (m)')
ax.set_title('3D Orbit of the Satellite')

# Zoom out
ax.set_xlim(-8000000, 8000000)
ax.set_ylim(-8000000, 8000000)
ax.set_zlim(-8000000, 8000000)

plt.show()