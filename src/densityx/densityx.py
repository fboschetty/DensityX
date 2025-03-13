import numpy as np
import pandas as pd


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
	from thermo_props import ThermodynamicProperties

	td_props = ThermodynamicProperties()

	data = dataframe.copy()
	mole_proportion = data / td_props.molecular_weight
	mole_fraction = mole_proportion.div(mole_proportion[OXIDE_COLUMNS].sum(axis=1), axis=0)
	return mole_fraction


def Density(dataframe: pd.DataFrame, verbose: bool = False) -> pd.DataFrame:
	"""Calculate density using model of Lange and Carmicheal 1990.

	Args:
		dataframe (pd.DataFrame): dataframe containing compositions in wt% oxides, Pressure in bar and temperature in celsius.
		verbose (bool, optional): flag, if True returns density and intermediate steps. Defaults to False.

	Raises:
		ValueError: if input is missing required columns.

	Returns:
		pd.DataFrame: calculated densities and uncertainty.
	"""
	from thermo_props import ThermodynamicProperties

	td_props = ThermodynamicProperties()

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
