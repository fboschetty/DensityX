import pytest
import numpy as np
import pandas as pd
from src.densityx import normalize_wt_percent_vals, mole_fraction, ThermodynamicProperties, Density


class TestThermodynamicProperties:
    @pytest.fixture
    def thermo_props(self):
        # Creates an instance of ThermodynamicProperties for use in tests
        return ThermodynamicProperties()

    def test_initialization(self, thermo_props):
        # Test that the class is initialized with the correct values

        # Test if 'molecular_weight' is a pandas Series
        assert isinstance(thermo_props.molecular_weight, pd.Series)

        # Test if 'molecular_weight' has the expected value for SiO2
        assert thermo_props.molecular_weight["SiO2"] == 60.0855

        # Test other attributes in a similar way
        assert thermo_props.molar_volume["TiO2"] == 28.32
        assert thermo_props.uncertainty_molar_volume["Al2O3"] == 0.09
        assert thermo_props.thermal_expansion_coefficient["H2O"] == 0.0095
        assert thermo_props.compressibility["Na2O"] == -0.00024
        assert thermo_props.reference_temperature["CaO"] == 1773

    def test_undefined_values(self, thermo_props):
        # Test if undefined properties are handled as expected (0.0 or 0)

        # For example, SiO2 has no thermal expansion coefficient, so we expect it to be 0.0
        assert thermo_props.thermal_expansion_coefficient["SiO2"] == 0.0

        # Similarly, check for other attributes like uncertainty values
        assert thermo_props.uncertainty_thermal_expansion["Fe2O3"] == 0
        assert thermo_props.uncertainty_compressibility["MgO"] == 0.000007

    def test_invalid_oxide(self, thermo_props):
        # Test what happens when an invalid oxide is queried
        with pytest.raises(KeyError):
            thermo_props.molecular_weight["InvalidOxide"]

        with pytest.raises(KeyError):
            thermo_props.molar_volume["InvalidOxide"]

    def test_property_values(self, thermo_props):
        # Test that the values of thermodynamic properties are within a reasonable range

        # Check if molecular weights are positive
        for value in thermo_props.molecular_weight:
            assert value > 0

        # Check if molar volumes are positive (e.g., in cm³/mol)
        for value in thermo_props.molar_volume:
            assert value > 0

        # Check if uncertainties are non-negative
        for value in thermo_props.uncertainty_molar_volume:
            assert value >= 0

        # Check if thermal expansion coefficients are reasonable (could be positive or zero)
        for value in thermo_props.thermal_expansion_coefficient:
            assert value >= 0.0

    def test_reference_temperature(self, thermo_props):
        # Test the reference temperature values
        assert thermo_props.reference_temperature["H2O"] == 1273
        assert thermo_props.reference_temperature["TiO2"] == 1773


class TestNormalizeWtPercentVals:
    def test_normalization_basic(self):
        data = pd.DataFrame({
            "A": [50, 30],
            "B": [50, 70]
        })
        result = normalize_wt_percent_vals(data)
        expected = pd.DataFrame({
            "A": [50.0, 30.0],
            "B": [50.0, 70.0],
            "Sum": [100.0, 100.0]
        })
        pd.testing.assert_frame_equal(result, expected)


    def test_normalization_non_100_totals(self):
        data = pd.DataFrame({
            "A": [25, 60],
            "B": [75, 40]
        })
        result = normalize_wt_percent_vals(data)
        expected = pd.DataFrame({
            "A": [25.0, 60.0],
            "B": [75.0, 40.0],
            "Sum": [100.0, 100.0]
        })
        pd.testing.assert_frame_equal(result, expected)


    def test_normalization_with_zeroes(self):
        data = pd.DataFrame({
            "A": [0., 50.],
            "B": [0., 50.]
        })
        result = normalize_wt_percent_vals(data)
        expected = pd.DataFrame({
            "A": [np.nan, 50.0],
            "B": [np.nan, 50.0],
            "Sum": [0.0, 100.0]  # First row remains 0 because sum was 0
        })
        pd.testing.assert_frame_equal(result, expected)


    def test_normalization_with_nan(self):
        data = pd.DataFrame({
            "A": [50, None],
            "B": [50, 100]
        })
        result = normalize_wt_percent_vals(data)
        expected = pd.DataFrame({
            "A": [50.0, None],
            "B": [50.0, 100.0],
            "Sum": [100.0, 100.0]
        })
        pd.testing.assert_frame_equal(result, expected)


class TestMoleFraction:
    @pytest.fixture
    def td_props(self):
        # Create an instance of ThermodynamicProperties for use in tests
        return ThermodynamicProperties()

    @pytest.fixture
    def sample_dataframe(self):
        # Sample dataframe to be used in tests
        data = {
            "SiO2": [60, 30],
            "TiO2": [15, 20],
            "Al2O3": [25, 50],
        }
        return pd.DataFrame(data)

    def test_mole_fraction_correctness(self, sample_dataframe, td_props):
        # Test that the mole_fraction function correctly computes mole fractions
        result = mole_fraction(sample_dataframe)

        # Calculate expected mole fractions for each oxide
        mole_fractions_SiO2 = sample_dataframe["SiO2"] / td_props.molecular_weight["SiO2"]
        mole_fractions_TiO2 = sample_dataframe["TiO2"] / td_props.molecular_weight["TiO2"]
        mole_fractions_Al2O3 = sample_dataframe["Al2O3"] / td_props.molecular_weight["Al2O3"]

        # Sum mole fractions for all oxides in each row
        total_mole_fractions = mole_fractions_SiO2 + mole_fractions_TiO2 + mole_fractions_Al2O3

        # Normalize each oxide's mole fraction by the total mole fraction for each row
        expected_mole_fraction_SiO2 = mole_fractions_SiO2 / total_mole_fractions
        expected_mole_fraction_TiO2 = mole_fractions_TiO2 / total_mole_fractions
        expected_mole_fraction_Al2O3 = mole_fractions_Al2O3 / total_mole_fractions

        # Assert that the results match the expected mole fractions
        assert np.allclose(result["SiO2"], expected_mole_fraction_SiO2)
        assert np.allclose(result["TiO2"], expected_mole_fraction_TiO2)
        assert np.allclose(result["Al2O3"], expected_mole_fraction_Al2O3)

    def test_empty_dataframe(self):
        # Test if mole_fraction handles an empty dataframe
        empty_df = pd.DataFrame()
        result = mole_fraction(empty_df)
        assert result.empty  # The result should be an empty dataframe

    def test_missing_columns(self, sample_dataframe):
        # Test for missing necessary oxide columns
        missing_columns_df = sample_dataframe.drop(columns=["TiO2"])
        with pytest.raises(KeyError):
            mole_fraction(missing_columns_df)

    def test_shape_consistency(self, sample_dataframe):
        # Ensure that the shape of the returned dataframe is the same as the input
        result = mole_fraction(sample_dataframe)
        assert result.shape == sample_dataframe.shape

    def test_nan_and_zero_values(self, sample_dataframe, td_props):
        # Test that NaN or zero values in the input dataframe are handled correctly
        sample_dataframe_with_nan = sample_dataframe.copy()
        sample_dataframe_with_nan.loc[0, "SiO2"] = np.nan
        sample_dataframe_with_nan.loc[1, "Al2O3"] = 0

        result = mole_fraction(sample_dataframe_with_nan)

        # Check that NaN or zero values are handled appropriately
        assert np.isnan(result["SiO2"].iloc[0])  # NaN should propagate in mole fractions
        assert result["Al2O3"].iloc[1] == 0  # Zero should result in a mole fraction of 0


class TestDensity:
    @pytest.fixture
    def td_props(self):
        # Create an instance of ThermodynamicProperties for use in tests
        return ThermodynamicProperties()

    @pytest.fixture
    def sample_dataframe(self):
        # Sample dataframe to be used in tests
        data = {
            "Sample_ID": [1, 2],
            "SiO2": [40.30, 45.60],
            "TiO2": [0.56, 3.63],
            "Al2O3": [16.24, 15/97],
            "Fe2O3": [0., 0.],
            "FeO": [5.59, 9.48],
            "MgO": [0.73, 4.28],
            "CaO": [18.74, 15.68],
            "Na2O": [12.23, 3.32],
            "K2O": [5.22, 0.32],
            "H2O": [0., 0.,],
            "P": [4000, 5000],  # Pressure in bar
            "T": [1000, 1200],  # Temperature in celsius
        }
        return pd.DataFrame(data)

    def test_missing_columns(self, sample_dataframe):
        # Test if the function raises ValueError when required columns are missing
        missing_columns_df = sample_dataframe.drop(columns=["SiO2"])
        with pytest.raises(ValueError):
            Density(missing_columns_df)

    def test_density_calculation(self, sample_dataframe, td_props):
        # Test that the function calculates density correctly
        result = Density(sample_dataframe)

        # Check that the output contains the expected columns
        assert "density_g_per_cm3" in result.columns
        assert "density_unc_g_per_cm3" in result.columns

        # You can also verify against manually computed expected values here,
        # but given the complexity of the calculation, we may not test for exact values in this test.
        # Example assertion (to verify presence of calculated density):
        assert result["density_g_per_cm3"].iloc[0] > 0  # Ensure density is a positive value

    def test_density_verbose(self, sample_dataframe):
        # Test that the function returns verbose output when verbose=True
        result = Density(sample_dataframe, verbose=True)

        # Check that the verbose output contains all the intermediate steps
        assert "Norm_SiO2" in result.columns  # Normalized values should have a prefix 'Norm_'
        assert "MoleFrac_SiO2" in result.columns  # Mole fraction values should have a prefix 'MoleFrac_'
        assert "Compdensity_SiO2" in result.columns  # Component density values should have a prefix 'Compdensity_'
        assert "VLiq_SiO2" in result.columns  # Liquid molar volume values should have a prefix 'VLiq_'
        assert "XMW_SiO2" in result.columns  # X * MW values should have a prefix 'XMW_'
        assert "AbsError_SiO2" in result.columns  # Absolute error values should have a prefix 'AbsError_'

    def test_empty_dataframe(self):
        # Test if the function handles an empty dataframe
        empty_df = pd.DataFrame(columns=["Sample_ID", "SiO2", "TiO2", "Al2O3", "P", "T"])
        with pytest.raises(ValueError):
            Density(empty_df)

    def test_nan_and_zero_values(self, sample_dataframe):
        # Test how the function handles NaN or zero values
        sample_dataframe_with_nan = sample_dataframe.copy()
        sample_dataframe_with_nan.loc[0, "SiO2"] = np.nan
        sample_dataframe_with_nan.loc[1, "Al2O3"] = 0

        result = Density(sample_dataframe_with_nan)

        # Ensure that NaN and zero values are handled properly (e.g., density should not be NaN or zero)
        assert not np.isnan(result["density_g_per_cm3"].iloc[0])  # Density should not be NaN
        assert result["density_g_per_cm3"].iloc[1] > 0  # Ensure positive density even with zero values
        assert result["density_unc_g_per_cm3"].iloc[1] >= 0  # Ensure uncertainty is not negative

    def test_correct_columns_in_result(self, sample_dataframe):
        # Test that the required columns are present in the result
        result = Density(sample_dataframe)

        # Check for necessary columns in the result dataframe
        assert "Sample_ID" in result.columns
        assert "density_g_per_cm3" in result.columns
        assert "density_unc_g_per_cm3" in result.columns

    def test_large_input_dataframe(self):
        # Test if the function handles large input data gracefully
        large_data = {
            "Sample_ID": range(1000),
            "SiO2": np.random.rand(1000) * 100,
            "TiO2": np.random.rand(1000) * 100,
            "Al2O3": np.random.rand(1000) * 100,
            "Fe2O3": np.random.rand(1000) * 100,
            "FeO": np.random.rand(1000) * 100,
            "MgO": np.random.rand(1000) * 100,
            "CaO": np.random.rand(1000) * 100,
            "Na2O": np.random.rand(1000) * 100,
            "K2O": np.random.rand(1000) * 100,
            "H2O": np.random.rand(1000) * 100,
            "P": np.random.rand(1000) * 100,
            "T": np.random.rand(1000) * 1000,
        }
        large_df = pd.DataFrame(large_data)
        result = Density(large_df)

        # Ensure that the result is not empty and has the correct number of rows
        assert not result.empty
        assert result.shape[0] == 1000  # Same number of rows as input
