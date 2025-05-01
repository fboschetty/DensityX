import numpy as np
import pandas as pd
from densityx import Density, normalize_wt_percent_vals, mole_fraction, ThermodynamicProperties

s_data = {
        "Sample_ID": 1,
        "SiO2":  50.08,
        "TiO2":   1.84,
        "Al2O3": 13.70,
        "Fe2O3":  2.70,
        "FeO":    9.57,
        "MgO":    6.67,
        "CaO":   11.50,
        "Na2O":   2.68,
        "K2O":    0.25,
        "H2O":    0.00,
        "P":    500.00,  # Pressure in bar
        "T":   1200.00,  # Temperature in celsius
    }

dataframe = pd.DataFrame(s_data, index=[0])



OXIDE_COLUMNS = ["SiO2", "TiO2", "Al2O3", "Fe2O3", "FeO", "MgO", "CaO", "Na2O", "K2O", "H2O"]

td_props = ThermodynamicProperties()

missing_columns = [
    col for col in OXIDE_COLUMNS + ["Sample_ID", "P", "T"]
    if col not in dataframe.columns
    ]
if missing_columns:
    raise ValueError(f"The following columns are missing from the input: {missing_columns}")

data = dataframe.copy()
data = data.fillna(value=0)
data_oxides = data[OXIDE_COLUMNS]

normalized = normalize_wt_percent_vals(data_oxides)
mole_fraction_vals = mole_fraction(normalized[OXIDE_COLUMNS])

# Convert temperatures to Kelvin
# Ensure P and T are in mole_fraction for subsequent calculations
mole_fraction_vals["T_K"] = data["T"] + 273.15
mole_fraction_vals["P"] = data["P"]

numerator = mole_fraction_vals * td_props.molecular_weight
denominator = mole_fraction_vals.apply(
    lambda row: td_props.molar_volume + td_props.thermal_expansion_coefficient * (row["T_K"] - td_props.reference_temperature) + td_props.compressibility * (row["P"] - 1),
    axis=1
    )
component_density = mole_fraction_vals[OXIDE_COLUMNS] * td_props.molecular_weight / denominator[OXIDE_COLUMNS]

# Calculate liquid molar volumes
Vliq = mole_fraction_vals * mole_fraction_vals.apply(
    lambda row: td_props.molar_volume + td_props.thermal_expansion_coefficient * (row["T_K"] - td_props.reference_temperature) + td_props.compressibility * (row["P"] - 1),
    axis=1
    )
Vliq["Sum"] = Vliq[OXIDE_COLUMNS].sum(axis=1)

# Calculate X*MW
X_MW = mole_fraction_vals[OXIDE_COLUMNS] * td_props.molecular_weight
X_MW["Sum"] = X_MW[OXIDE_COLUMNS].sum(axis=1)

# Calculate the density of the melt in g/cm3 and in g/L
density_g_per_cm3 = X_MW["Sum"] / Vliq["Sum"]

print(f"Norm: {normalized[OXIDE_COLUMNS]}")
print(f"Mol_Frac: {mole_fraction_vals[OXIDE_COLUMNS]}")
print(f"Num: {numerator[OXIDE_COLUMNS]}")
print(f"Denom: {denominator[OXIDE_COLUMNS]}")
print(Vliq[OXIDE_COLUMNS])
print(Vliq["Sum"])

print(density_g_per_cm3)