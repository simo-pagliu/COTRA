import numpy as np
from scipy.integrate import solve_ivp
import os
import glob
import re

def run(Hydro_Dispersion, pore_velocity, porosity, bulk_density,
        Source_Time, Source_Intensity, Domain_Length, dx, Time_Span, dt,
        retardation, k_1, k_2):
    
    Grid_Space = np.arange(0, Domain_Length + dx, dx)
    Steps_Space = len(Grid_Space)
    Grid_Time = np.arange(Time_Span[0], Time_Span[1] + dt, dt)
    Concentration_Init = np.zeros(Steps_Space)
    Concentration_Init[0] = Source_Intensity

    def pde_rhs(t, C, source_active=True):
        C_ext = np.zeros(Steps_Space + 2)
        C_ext[1:-1] = C
        C_ext[0] = C[1]
        C_ext[-1] = C[-2]

        dC_dx = (C_ext[2:] - C_ext[:-2]) / (2 * dx)
        d2C_dx2 = (C_ext[2:] - 2 * C_ext[1:-1] + C_ext[:-2]) / (dx ** 2)

        if retardation == 1:
            R = 1 + k_1 * bulk_density / porosity
        elif retardation == 2:
            R = 1 + k_1 * bulk_density / porosity * (C ** (k_2 - 1))
        else:
            R = 1.0

        dC_dt = (-pore_velocity * dC_dx + Hydro_Dispersion * d2C_dx2) / R

        if source_active:
            C[0] = Source_Intensity
        else:
            C[0] = 0

        dC_dt[0] = 0
        dC_dt[-1] = 0
        return dC_dt

    t_span1 = (Time_Span[0], Source_Time)
    t_eval1 = Grid_Time[Grid_Time <= Source_Time]
    sol1 = solve_ivp(lambda t, C: pde_rhs(t, C, True), t_span1, Concentration_Init, t_eval=t_eval1, method='Radau')

    t_span2 = (Source_Time, Time_Span[1])
    t_eval2 = Grid_Time[Grid_Time >= Source_Time]
    Concentration_at_switch = sol1.y[:, -1].copy()
    sol2 = solve_ivp(lambda t, C: pde_rhs(t, C, False), t_span2, Concentration_at_switch, t_eval=t_eval2, method='RK45')

    t_full = np.concatenate((sol1.t, sol2.t))
    C_full = np.hstack((sol1.y, sol2.y))

    base_name, file_ext = "output_data", ".npz"
    pattern = base_name + "*" + file_ext
    existing_files = glob.glob(pattern)

    max_counter = -1
    for file in existing_files:
        match = re.match(rf"{base_name}(?:_(\d+))?{re.escape(file_ext)}", os.path.basename(file))
        if match:
            num_str = match.group(1)
            max_counter = max(max_counter, int(num_str) if num_str else 0)

    new_counter = max_counter + 1 if max_counter >= 0 else 0
    filename = f"{base_name}{file_ext}" if new_counter == 0 else f"{base_name}_{new_counter}{file_ext}"

    np.savez(filename, t_full=t_full, C_full=C_full, Domain_Length=Domain_Length,
             Source_Intensity=Source_Intensity, Grid_Space=Grid_Space)

    print(f"Data saved to {filename}")
    return filename
