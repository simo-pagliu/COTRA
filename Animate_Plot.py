import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Load the file
data = np.load("output_data.npz")

# Extract variables (adjust key names if needed)
t_full = data["t_full"]
C_full = data["C_full"]
Domain_Length = data["Domain_Length"].item()  # Convert scalar array to a Python scalar
Source_Intensity = data["Source_Intensity"].item()
Grid_Space = data["Grid_Space"]

# Animate the concentration profile evolution
fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
ax.set_xlim(0, Domain_Length)
ax.set_ylim(0, Source_Intensity * 1.1)
ax.set_xlabel('x [m]')
ax.set_ylabel('Concentration [g/l]')
ax.set_title('Concentration Profile Evolution')

def init():
    line.set_data([], [])
    return line,

def update(frame):
    line.set_data(Grid_Space, C_full[:, frame])
    ax.set_title(f'Concentration Profile at t = {t_full[frame]:.0f} s')
    return line,

frames = range(0, len(t_full), 10)  # Use every 10th time step to reduce frames
ani = FuncAnimation(fig, update, frames=frames, init_func=init, blit=True)

plt.show()