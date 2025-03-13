"""
DensityX

A python library for calculating the densities of silicate melts.

Published as: Iacovino, K. and Till, C. B. (2019) “DensityX: A program for calculating the
densities of magmatic liquids up to 1,627 °C and 30 kbar”, Volcanica, 2(1), pp. 1-10.
doi: 10.30909/vol.02.01.0110.
"""

__version__ = "1.2.0"
__author__ = "Kayla Iacovino, Christy Till, Felix Boschetty"


import numpy as np
import pandas as pd
from dataclasses import dataclass, field


class ThermodynamicProperties:
    """
    A class that holds the thermodynamic properties for various oxides.
	Undefined values are set to 0.

    Attributes:
        molecular_weight (pd.Series): Molecular weights for oxides given in g/mol.
        molar_volume (pd.Series): Molar volumes for oxides in cm³/mol.
        uncertainty_molar_volume (pd.Series): Reported uncertainties in molar volumes for oxides in cm³/mol.
        thermal_expansion_coefficient (pd.Series): Thermal expansion coefficient for each oxide in cm³/mol/K.
        uncertainty_thermal_expansion (pd.Series): Reported uncertainties in thermal expansion coefficient in cm³/mol/K.
        compressibility (pd.Series): Compressibility for each oxide given in cm³/mol/bar.
        uncertainty_compressibility (pd.Series): Reported uncertainties in compressibility for each oxide in cm³/mol/bar.
        reference_temperature (pd.Series): Reference temperatures in Kelvin for the above thermodynamic parameters.
    """

    def __init__(self):
        self.molecular_weight = pd.Series({
            "SiO2": 60.0855, "TiO2": 79.88, "Al2O3": 101.96, "Fe2O3": 159.69,
            "FeO": 71.85, "MgO": 40.3, "CaO": 56.08, "Na2O": 61.98, "K2O": 94.2, "H2O": 18.02,
        })

        self.molar_volume = pd.Series({
            "SiO2": 26.86, "TiO2": 28.32, "Al2O3": 37.42, "Fe2O3": 41.50,
            "FeO": 12.68, "MgO": 12.02, "CaO": 16.90, "Na2O": 29.65, "K2O": 47.28, "H2O": 22.9,
        })
        """Molar volumes for oxides in cm³/mol.

			* SiO₂, Al₂O₃, MgO, CaO, Na₂O, K₂O at 1773 K (Lange, 1997; CMP).
			* H₂O at 1273 K (Ochs and Lange, 1999).
			* FeO at 1723 K (Guo et al., 2014).
			* Fe₂O₃ at 1723 K (Liu and Lange, 2006).
			* TiO₂ at 1773 K (Lange and Carmichael, 1987).
        """

        self.uncertainty_molar_volume = pd.Series({
            "SiO2": 0.03, "TiO2": 0, "Al2O3": 0.09, "Fe2O3": 0,
            "FeO": 0, "MgO": 0.07, "CaO": 0.06, "Na2O": 0.07, "K2O": 0.10, "H2O": 0.60,
        })
        """Reported uncertainties in molar volumes for oxides in cm³/mol. References same as `molar_volumes`."""

        self.thermal_expansion_coefficient = pd.Series({
            "SiO2": 0.0, "TiO2": 0.00724, "Al2O3": 0.00262, "Fe2O3": 0.0,
            "FeO": 0.00369, "MgO": 0.00327, "CaO": 0.00374, "Na2O": 0.00768, "K2O": 0.01208, "H2O": 0.0095,
        })
        """Thermal expansion coefficient for each oxide in cm³/mol/K."""

        self.uncertainty_thermal_expansion = pd.Series({
            "SiO2": 0, "TiO2": 0, "Al2O3": 0, "Fe2O3": 0,
            "FeO": 0, "MgO": 0, "CaO": 0, "Na2O": 0, "K2O": 0, "H2O": 0.00080,
        })
        """Reported uncertainties in thermal expansion coefficient in cm³/mol/K."""

        self.compressibility = pd.Series({
            "SiO2": -0.000189, "TiO2": -0.000231, "Al2O3": -0.000226, "Fe2O3": -0.000253,
            "FeO": -0.000045, "MgO": 0.000027, "CaO": 0.000034, "Na2O": -0.00024,
            "K2O": -0.000675, "H2O": -0.00032,
        })
        """Compressibility for each oxide given in cm³/mol/bar."""

        self.uncertainty_compressibility = pd.Series({
            "SiO2": 0.000002, "TiO2": 0.000006, "Al2O3": 0.000009, "Fe2O3": 0.000009,
            "FeO": 0.000003, "MgO": 0.000007, "CaO": 0.000005, "Na2O": 0.000005,
            "K2O": 0.000014, "H2O": 0.000060,
        })
        """Reported uncertainties in compressibility for each oxide in cm³/mol/bar."""

        self.reference_temperature = pd.Series({
            "SiO2": 1773, "TiO2": 1773, "Al2O3": 1773, "Fe2O3": 1723,
            "FeO": 1723, "MgO": 1773, "CaO": 1773, "Na2O": 1773, "K2O": 1773, "H2O": 1273,
        })
        """Reference temperatures in Kelvin for thermodynamic parameters."""

td_props = ThermodynamicProperties()

OXIDE_COLUMNS = ["SiO2", "TiO2", "Al2O3", "Fe2O3", "FeO", "MgO", "CaO", "Na2O", "K2O", "H2O"]

def normalize_wt_percent_vals(dataframe: pd.DataFrame) -> pd.DataFrame:
	"""Normalize a dataframe of input weight percent values by 100 wt.% total.

	Args:
		dataframe (pd.DataFrame): input dataframe.

	Returns:
		pd.DataFrame: normalized dataframe.
	"""
	data = dataframe.copy()
	total = data.sum(axis=1)
	normalized = 100. * data.div(total, axis=0)
	normalized["Sum"] = normalized.sum(axis=1)
	return normalized


def mole_fraction(dataframe: pd.DataFrame) -> pd.DataFrame:
	"""Calculate the mole fractions of an input dataframe in oxide wt%.

	Args:
		dataframe (pd.DataFrame): input dataframe in oxide wt.%.

	Returns:
		pd.DataFrame: mole fractions.
	"""
	data = dataframe.copy()
	mole_proportion = data / td_props.molecular_weight
	mole_fraction = mole_proportion.div(mole_proportion[OXIDE_COLUMNS].sum(axis=1), axis=0)
	return mole_fraction


def density(dataframe: pd.DataFrame, verbose: bool = False) -> pd.DataFrame:
	"""Calculate density using model of Lange and Carmicheal 1990.

	Args:
		dataframe (pd.DataFrame): dataframe containing compositions in wt% oxides, Pressure in bar and temperature in celsius.
		verbose (bool, optional): flag, if True returns density and intermediate steps. Defaults to False.

	Raises:
		ValueError: if input is missing required columns.

	Returns:
		pd.DataFrame: calculated densities and uncertainty.
	"""
	missing_columns = [
		col for col in OXIDE_COLUMNS+["Sample_ID", "P", "T"]
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

	# Calculate component density.
	numerator = mole_fraction_vals * td_props.molecular_weight
	denominator = mole_fraction_vals.apply(
		lambda row: td_props.molecular_weight + td_props.thermal_expansion_coefficient * (row["T_K"] - td_props.reference_temperature) + td_props.compressibility * (row["P"] - 1),
		axis=1
		)
	component_density = numerator[OXIDE_COLUMNS] / denominator[OXIDE_COLUMNS]

	# Calculate liquid molar volumes
	Vliq = component_density * mole_fraction_vals.apply(
		lambda row: td_props.molar_volume + td_props.thermal_expansion_coefficient * (row["T_K"] - td_props.reference_temperature) + td_props.compressibility * (row["P"] - 1),
		axis=1
		)
	Vliq["Sum"] = Vliq[OXIDE_COLUMNS].sum(axis=1)

	# Calculate X*MW
	X_MW = mole_fraction_vals[OXIDE_COLUMNS] * td_props.molecular_weight
	X_MW["Sum"] = X_MW[OXIDE_COLUMNS].sum(axis=1)

	# Calculate the density of the melt in g/cm3 and in g/L
	density_g_per_cm3 = X_MW["Sum"] / Vliq["Sum"]

	# Calculate relative errors -> replace divide by zero errors by 0.
	rel_error_mv = (td_props.uncertainty_molar_volume / td_props.molar_volume).fillna(value=0.)
	rel_error_te = (td_props.uncertainty_thermal_expansion / td_props.thermal_expansion_coefficient).fillna(value=0.)
	rel_error_compress = (td_props.uncertainty_compressibility / td_props.compressibility).fillna(value=0.)

	relative_error_Vliq = np.sqrt(rel_error_mv**2 + rel_error_te**2 + rel_error_compress**2)
	absolute_error_Vliq = relative_error_Vliq * Vliq
	absolute_error_Vliq["Sum"] = absolute_error_Vliq[OXIDE_COLUMNS].sum(axis=1)

	# Calculate error on density value
	uncertainty_g_per_cm3 = absolute_error_Vliq["Sum"] / Vliq["Sum"]

	if verbose:
		data_to_return = pd.concat([
			data,
			normalized.add_prefix(prefix="Norm_"),
			mole_fraction_vals[OXIDE_COLUMNS].add_prefix(prefix="MoleFrac_"),
			component_density.add_prefix(prefix="Compdensity_"),
			Vliq.add_prefix(prefix="VLiq_"),
			X_MW.add_prefix(prefix="XMW_"),
			absolute_error_Vliq.add_prefix(prefix="AbsError_")
			], axis=1
		)
		data_to_return["density_g_per_cm3"] = density_g_per_cm3
		data_to_return["density_unc_g_per_cm3"] = uncertainty_g_per_cm3
		data_to_return["density_g_per_L"] = density_g_per_cm3 * 1000
		data_to_return["uncertainty_g_per_L"] = uncertainty_g_per_cm3 * 1000

	else:
		data_to_return = pd.DataFrame(data={
		"Sample_ID": data["Sample_ID"],
		"density_g_per_cm3": density_g_per_cm3,
		"density_unc_g_per_cm3": uncertainty_g_per_cm3
		})

	return data_to_return
