import pandas as pd


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