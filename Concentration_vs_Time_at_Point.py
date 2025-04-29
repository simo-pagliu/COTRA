import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# List of files to compare and their comments
# files = ["output_data_1.npz", "output_data_2.npz", "output_data_3.npz", "output_data_4.npz"]
# comments = ["No Retardation + x10 Hyd.Dispersion",
#             "No Retardation",
#             "Linear Retardation",
#             "Linear Retardation + x10 Hyd.Dispersion"]

files = ["output_data.npz"]
comments = ["simulation"]

# List to hold the loaded data
data_list = []

# Load each file and extract variables
for f in files:
    data = np.load(f)
    t_full = data["t_full"]
    C_full = data["C_full"]
    Domain_Length = data["Domain_Length"].item()  # Convert scalar array to a Python scalar
    Source_Intensity = data["Source_Intensity"].item()
    Grid_Space = data["Grid_Space"]
    data_list.append((t_full, C_full, Domain_Length, Source_Intensity, Grid_Space))

# Get array of concentration in time at a given point in space
point = 0.29  # Example point in space (in meters)

# Assume all files share the same time and space grids
t_full = data_list[0][0]
Grid_Space = data_list[0][4]
num_time_steps = len(t_full)
min_x = Grid_Space.min()
max_x = Grid_Space.max()

# Create the figure and the main axes for plotting
fig, ax = plt.subplots(figsize=(8, 5))
plt.subplots_adjust(bottom=0.25)  # Make space at the bottom for the slider

# Plot concentration in time at the selected spatial point for each dataset
lines = []
for idx, (t_full, C_full, Domain_Length, Source_Intensity, Grid_Space) in enumerate(data_list):
    point_index = np.argmin(np.abs(Grid_Space - point))
    concentration_evolution = C_full[point_index, :]  # Concentration over time at this point
    line, = ax.plot(t_full, concentration_evolution, label=f'{comments[idx]}')
    lines.append(line)

ax.set_xlabel('Time [s]')
ax.set_ylabel('Concentration [g/l]')
ax.set_title(f'Concentration vs Time at x = {point:.2f} m')
ax.legend()

# Create a slider axis below the main plot
slider_ax = plt.axes([0.15, 0.1, 0.7, 0.05])
x_slider = Slider(slider_ax, 'x [m]', min_x, max_x, valinit=point)

# Update function to change the plot based on the slider's value
def update(val):
    point = x_slider.val
    point_index = np.argmin(np.abs(Grid_Space - point))
    for idx, (_, C_full, _, _, _) in enumerate(data_list):
        concentration_evolution = C_full[point_index, :]
        lines[idx].set_ydata(concentration_evolution)
    ax.set_title(f'Concentration vs Time at x = {point:.2f} m')
    fig.canvas.draw_idle()

# Register the update function with the slider
x_slider.on_changed(update)

plt.show()

# Save the points of the profile at a given point in space
point = 19.5e-2
point_index = np.argmin(np.abs(Grid_Space - point))
concentration_evolution = C_full[point_index, :]
# Save as csv
np.savetxt("concentration_evolution.csv", concentration_evolution, delimiter=",")