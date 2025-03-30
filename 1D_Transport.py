import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.animation import FuncAnimation
import os
import glob
import re

# Input Parameters
Hydro_Dispersion = 3.3e-7       # Hydrodynamic dispersion [m^2/s]
pore_velocity = 1.0e-4          # Pore velocity [m/s]
porosity = 0.30                 # Porosity
bulk_density = 4                # Bulk density [kg/l]

# Source term parameters
Source_Time = 10000             # Source duration [s]
Source_Intensity = 1            # Source intensity [g/l]

# Discretization Settings
Domain_Length = 1               # Domain length [m]
dx = 0.005                    # Space discretization [m]
Time_Span = (0, 20000)          # Total simulation time [s]
dt = 40                       # Time step for evaluation [s]

# Create spatial grid
Grid_Space = np.arange(0, Domain_Length + dx, dx)
Steps_Space = len(Grid_Space)

# Create temporal grid for evaluation (for the entire simulation)
Grid_Time = np.arange(Time_Span[0], Time_Span[1] + dt, dt)

# Initial condition: a step function at the left boundary
Concentration_Init = np.zeros(Steps_Space)
Concentration_Init[0] = Source_Intensity

# Retardation settings (example, adjust as needed)
retardation = 0  # Options: 0, 1, or 2
k_1 = 0.1
k_2 = 0.2

def pde_rhs(t, C, source_active=True):
    """
    Right-hand side for the transport PDE using the method-of-lines.
    PDE:
        Retardation_Factor * dC/dt = pore_velocity * dC/dx + Hydro_Dispersion * d2C/dx2
    Forcing is applied at the left boundary when source_active is True.
    """
    # Compute spatial derivatives using finite differences
    dC_dx = np.gradient(C, dx)
    d2C_dx2 = np.gradient(dC_dx, dx)
    
    # Compute the retardation factor (here a constant or dependent on C)
    if retardation == 1:
        R = 1 + k_1 * bulk_density / porosity  # linear model
    elif retardation == 2:
        R = 1 + k_1 * bulk_density / porosity * (C ** (k_2 - 1))  # Freundlich model (use ** for exponentiation)
    else:
        R = 1.0  # no retardation

    # Compute the right-hand side of the PDE
    dC_dt = (-pore_velocity * dC_dx + Hydro_Dispersion * d2C_dx2) / R

    # If the source is active, enforce the boundary condition at x=0.
    # Here we simply force the derivative at the first node to be zero so that the value remains constant.
    if source_active:
        C[0] = Source_Intensity  # This holds C[0] constant
    else:
        C[0] = 0 
    return dC_dt

# First stage: with source active (t from 0 to Source_Time)
t_span1 = (Time_Span[0], Source_Time)
t_eval1 = Grid_Time[Grid_Time <= Source_Time]

sol1 = solve_ivp(lambda t, C: pde_rhs(t, C, source_active=True),
                 t_span1, Concentration_Init, t_eval=t_eval1, method='RK45')

# The final state from the first stage is used as the initial condition for the second stage.
Concentration_at_switch = sol1.y[:, -1].copy()

# Second stage: with the source switched off (t from Source_Time to final time)
t_span2 = (Source_Time, Time_Span[1])
t_eval2 = Grid_Time[Grid_Time >= Source_Time]

sol2 = solve_ivp(lambda t, C: pde_rhs(t, C, source_active=False),
                 t_span2, Concentration_at_switch, t_eval=t_eval2, method='RK45')

# Combine the time and solution arrays from both stages
t_full = np.concatenate((sol1.t, sol2.t))
C_full = np.hstack((sol1.y, sol2.y))

#################################################################
# Save the data to a file
#################################################################
base_name = "output_data"
file_ext = ".npz"
pattern = base_name + "*" + file_ext

# Find all files matching the pattern
existing_files = glob.glob(pattern)

max_counter = -1
for file in existing_files:
    # Expected filenames:
    # "output_data.npz" or "output_data_#.npz"
    base_file = os.path.basename(file)
    match = re.match(rf"{base_name}(?:_(\d+))?{re.escape(file_ext)}", base_file)
    if match:
        num_str = match.group(1)
        if num_str is None:
            # "output_data.npz" exists; treat its counter as 0.
            max_counter = max(max_counter, 0)
        else:
            max_counter = max(max_counter, int(num_str))

# Determine the new counter
new_counter = max_counter + 1 if max_counter >= 0 else 0

# Build the filename based on the counter.
if new_counter == 0:
    filename = base_name + file_ext
else:
    filename = f"{base_name}_{new_counter}{file_ext}"

# Save the data into the file.
np.savez(filename, t_full=t_full, C_full=C_full,
         Domain_Length=Domain_Length, Source_Intensity=Source_Intensity,
         Grid_Space=Grid_Space)

print(f"Data saved to {filename}")