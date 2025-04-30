# Importing the COTRA module
import COTRA

###############################################################################
# Parameters needed to solve the equation
# Matrix properties
Hydro_Dispersion = 1.6306*1e-6       # Hydrodynamic dispersion [m^2/s]
pore_velocity = 0.9734*1e-3          # Pore velocity [m/s]
porosity = 0.4078                 # Porosity
bulk_density = 1.48               # Bulk density [kg/l]

# Source term parameters
Source_Time = 480             # Source duration [s]
Source_Intensity = 1            # Source intensity [g/l]

# Discretization Settings
Domain_Length = 19.5e-2               # Domain length [m]
dx = 1e-3                    # Space discretization [m]
Time_Span = (0, 1200)          # Total simulation time [s]
dt = 1                       # Time step for evaluation [s]

# Retardation settings
# 0: no retardation, 
# 1: linear retardation,
# 2: nonlinear retardation (Freudlich model)
retardation = 0 # 0,1 or 2
k_1 = 0.1
k_2 = 0.2
###############################################################################

# Runs the code to solve the advection-dispersion equation
# The output is saved in "output_data.npz" 
# If the "output_data.npz" already exists it creates "output_data_1.npz", "output_data_2.npz", etc.
COTRA.run(Hydro_Dispersion, pore_velocity, porosity, bulk_density,
        Source_Time, Source_Intensity, Domain_Length, dx, Time_Span, dt,
        retardation, k_1, k_2)

# Plot of the concentartion vs space animated in time
COTRA.animate("output_data.npz")

# Interactive plot of concentration against space at different points in time
# The legend will be populated with "Example of COTRA.CvS"
COTRA.CvS("output_data.npz", "Example of COTRA.CvS")

# Interactive plot of concentration against time at different points in space
# No comments added: legend will be empty
COTRA.CvT("output_data.npz")

# Save concentration profile vs Time at a specific space point: 19.5 cm
COTRA.save_profile("output_data.npz", 19.5e-2)

# Save concentration profile vs Space at a specific time point: 150 s
COTRA.save_profile("output_data.npz", 150, domain='time')