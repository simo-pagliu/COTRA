1: No Retardation + x10 Hydrodynmaic Dispersion Coefficient
2: No Retardation
3: Linear Retardation
4: Linear Retardation + x10 Hydrodynmaic Dispersion Coefficient


Per il NON EQUILIBRIUM CASE (1 region)
cambia absorption rate ad un valore non nullo

da: dC_dt = (-pore_velocity * dC_dx + Hydro_Dispersion * d2C_dx2) / R

passiamo a:
dC_dt = (-pore_velocity * dC_dx + Hydro_Dispersion * d2C_dx2 - S) / R
density* ds_dt = adsorption_Rate * porosity * C - desorption_Rate * density * s

adsorption_rate = desorption_rate * density / porosity * k_1 * C ** (k_2 - 1) FREUDLICH
adsorption_rate = desorption_rate * density / porosity * k_1 LINEAR


TWO REGION MODEL (alternativa non equilibrium case)
MOBILE: dC_dt = (-pore_velocity * dC_dx + Hydro_Dispersion * d2C_dx2 - alpha * (C - C_immobile)) / R
IMMOBILE: dC_immobile_dt = (alpha * (C - C_immobile))/R

