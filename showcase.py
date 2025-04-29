import COTRA

# Input Parameters
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

# Retardation settings (example, adjust as needed)
retardation = 0  # Options: 0, 1, or 2
k_1 = 0.1
k_2 = 0.2


COTRA.run(Hydro_Dispersion, pore_velocity, porosity, bulk_density,
        Source_Time, Source_Intensity, Domain_Length, dx, Time_Span, dt,
        retardation, k_1, k_2)

COTRA.animate("output_data_5.npz")
COTRA.CvS("output_data_5.npz", "output_data_5.npz")
COTRA.CvT("output_data_5.npz", "output_data_5.npz")
COTRA.save_profile("output_data_5.npz", 19.5)