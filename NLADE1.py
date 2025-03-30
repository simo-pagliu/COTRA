import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# Input Parameters
Hydro_Dispersion = 3.3e-7       # Hydrodinamic dispersion [m^2/s]
pore_velocity = 1.0e-4       # pore velocity [m/s]
diffusion_rate_bw_regions = 0.0001    # diffusion rate between regions [1/s]
porosity = 0.30     # porosity
porosity_mobile = 0.25    # mobile porosity
porosity_immobile = 0.35   # immobile porosity
bulk_density = 4          # bulk density [kg/l]
# Source term parameters
Source_Time = 10000       # source duration [s]
Source_Intensity = 1           # source intensity [g/l]

# Discretization Settings
Domain_Length = 1            # domain length [m]
dx = 0.005       # space discretization [m]
Time_Span = (0, 20000)      # final simulation time [s]
dt = 40        # time step [s]

# Create spatial grid
Grid_Space = np.arange(0, Domain_Length + dx, dx)
Steps_Space = len(Grid_Space)  # number of spatial steps

# Creates temporal grid
Grid_Time = np.arange(Time_Span[0], Time_Span[1] + dt, dt)
Steps_Time = len(Grid_Time)  # number of time steps

# Initial conditions
source_stop_time_index = np.where(Grid_Time >= Source_Time)[0][0]
Concentration_Init = np.zeros((Steps_Space))
Concentration_Init[0] = Source_Intensity

retardation = 0 #1, 2
k_1 = 0.1
k_2 = 0.5
def pde_rhs(t, C):
    """
    Right-hand side for the transport PDE in method-of-lines form.
    
    The PDE is:
    
        Retardation_Factor * dC/dt = pore_velocity * dC/dx + Hydro_Dispersion * d2C/dx2
    
    Thus:
    
        dC/dt = [pore_velocity * dC/dx + Hydro_Dispersion * d2C/dx2] / Retardation_Factor
    
    Finite differences (using np.gradient) approximate the spatial derivatives.
    """
    # Compute the first derivative dC/dx (using central differences)
    dC_dx = np.gradient(C, dx)
    # Compute the second derivative d2C/dx2
    d2C_dx2 = np.gradient(dC_dx, dx)

    if retardation == 1:
        Retardation_Factor = 1 + k_1 * bulk_density / porosity # Linear model
    elif retardation == 2:
        Retardation_Factor = 1 + k_1 * bulk_density / porosity * C ^ (k_2 - 1) # Fraulich model
    else:
        Retardation_Factor = 1 # where 1 means no retardation
    
    # Compute time derivative dC/dt pointwise
    dC_dt = (pore_velocity * dC_dx + Hydro_Dispersion * d2C_dx2) / Retardation_Factor
    return dC_dt

# Solve the PDE (method-of-lines) using an ODE solver
print("Dimension of Concentration_Init:", Concentration_Init.shape)
sol = solve_ivp(pde_rhs, t_span = Time_Span, y0 = 
                Concentration_Init, t_eval=Grid_Time, method='RK45')

# Plot the concentration profile at the final time
plt.plot(Grid_Space, sol.y[:, -1], label=f't = {Grid_Time[-1]:.0f} s')
plt.xlabel('x')
plt.ylabel('Concentration')
plt.title('Concentration profile at final time')
plt.legend()
plt.show()



# # Find the index where the time grid first reaches or exceeds the source duration Source_Time
# # MATLAB: isend = min(find(Grid_Time>=Source_Time));
# 

# if par1 == 1 or par1 == 2:
#     # Simulation for linear (par1==1) or non-linear (par1==2) isotherm in 1 or 2 steps:
#     # First step: with source (time up to Source_Time)
#     T = Grid_Time[:isend+1]  # Include the time corresponding to Source_Time
#     out1 = pdeade(T, kdrt, n, Source_Intensity, Hydro_Dispersion, pore_velocity, C0nl, Grid_Space)
#     outi = out1.copy()

#     if isend < len(Grid_Time) - 1:
#         # Update initial condition from the last time step of out1
#         C0nl = np.squeeze(out1[-1, :]).reshape(-1, 1)
#         # Second step: simulation with source turned off (Source_Intensity = 0)
#         T = Grid_Time[isend:]
#         Source_Intensity = 0
#         out2 = pdeade(T, kdrt, n, Source_Intensity, Hydro_Dispersion, pore_velocity, C0nl, Grid_Space)
#         # Concatenate the results, skipping the duplicate first row of out2
#         outi = np.vstack([outi, out2[1:, :]])

#     # Save result to ASCII file (similar to MATLAB's "save ... -ascii")
#     np.savetxt('out_A.dat', outi, fmt='%g')

# elif par1 == 3 or par1 == 4:
#     # Simulation for sorption/double porosity cases (linear par1==3 or non-linear par1==4)
#     T = Grid_Time[:isend+1]
#     out1 = pdeade2(T, diffusion_rate_bw_regions, adsorption_rate, porosity, porosity_mobile, porosity_immobile, bulk_density, Source_Intensity, kdrt, n, Hydro_Dispersion, pore_velocity, C0nl, Grid_Space)
#     outi = out1.copy()

#     if isend < len(Grid_Time) - 1:
#         C0nl = np.squeeze(out1[-1, :, :])
#         T = Grid_Time[isend:]
#         Source_Intensity = 0
#         out2 = pdeade2(T, diffusion_rate_bw_regions, adsorption_rate, porosity, porosity_mobile, porosity_immobile, bulk_density, Source_Intensity, kdrt, n, Hydro_Dispersion, pore_velocity, C0nl, Grid_Space)
#         outi = np.concatenate((outi, out2[1:, :, :]), axis=0)

#     # Save results for the two particles (pages 1 and 2)
#     out_particle1 = outi[:, :, 0]
#     np.savetxt('out_1A.dat', out_particle1, fmt='%g')
#     out_particle2 = outi[:, :, 1]
#     np.savetxt('out_2A.dat', out_particle2, fmt='%g')

# # Plotting the results
# plt.figure()
# # Find the first spatial index where Grid_Space >= 0.1
# idx = np.where(Grid_Space >= 0.1)[0][0]
# if par1 == 1 or par1 == 2:
#     plt.plot(Grid_Time, outi[:, idx] / C0nl[0, 0])
# else:
#     plt.plot(Grid_Time, outi[:, idx, 0] / C0nl[0, 0])
# plt.ylabel('C/C0')
# plt.xlabel('t[s]')
# plt.show()

# print(outi, Grid_Space, Grid_Time)