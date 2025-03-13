"""
DensityX

A python library for calculating the densities of silicate melts.

Published as: Iacovino, K. and Till, C. B. (2019) “DensityX: A program for calculating the
densities of magmatic liquids up to 1,627 °C and 30 kbar”, Volcanica, 2(1), pp. 1-10.
doi: 10.30909/vol.02.01.0110.
"""

__version__ = "1.2.0"
__author__ = "Kayla Iacovino, Christy Till, Felix Boschetty"

from .densityx import normalize_wt_percent_vals, mole_fraction, Density
from .thermo_props import ThermodynamicProperties

# Define accessible objects when importing from the module
__all__ = [
    "ThermodynamicProperties",
    "normalize_wt_percent_vals",
    "mole_fraction",
    "Density"
]
