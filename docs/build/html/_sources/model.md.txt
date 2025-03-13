# DensityX Model

**DensityX** is a standalone Python package for calculating the density of silicate melts based on their major oxide compositions, pressure, and temperature. It relies on partial molar volumes of oxide components and assumes ideal mixing. The model utilises the equation of state of {cite:t}`lange1997` and is calibrated for a wide range of silicate melt compositions, operating over pressures of 1–30 kbar and temperatures up to 1627°C. DensityX does not require knowledge of the melt oxygen fugacity. Further details about the model are discussed below and in {cite:t}`iacovino2018`.

## Thermodynamic Data Sources

DensityX is built upon established thermodynamic data from:

* **{cite:t}`lange1987`:** Partial molar volumes of silicate melt components.
* **{cite:t}`lange1997`:** Temperature and pressure dependencies of silicate melt properties.
* **{cite:t}`kress1991`:** Redox effects and thermodynamic properties of Fe-bearing melts.
* **{cite:t}`ochs1999`:** Effect of H₂O on melt density.
* **{cite:t}`guo2014`:** Updated calibration of thermodynamic parameters.

## Mathematical Model

Density is calculated using the **ideal mixing assumption**:

```{math}
\rho_{\textrm{liq}}​=V_{\textrm{liq}}​\sum{X_i​ M_i​}​
```

where:

* {math}`\rho_{\textrm{liq}}`​ is density of the silicate liquid,
* {math}`X_i`​ is mole fraction of oxide component {math}`i`,
* {math}`M_i` is molecular weight of oxide component {math}`i`, and
* {math}`V_{\textrm{liq}}`​ is total volume of the silicate liquid at a given {math}`P` and {math}`T`, computed as:

```{math}
V_{\textrm{liq}} = \sum{X_i V_i}
```

where {math}`V_i`​ is the partial molar volume of each oxide component, adjusted for pressure and temperature using:

```{math}
V_{i}(T, P) = V_{i}(T_{\textrm{ref}}, P_{\textrm{ref}}) + dV_i dT(T−T_{\textrm{ref}}) + dV_{i}dP(P−P_{\textrm{ref}})
```

where {math}`T_{\textrm{ref}}` is the reference temperature for the molar volume, and {math}`P_{\textrm{ref}}` is the reference pressure, often 1 bar.

This equation is valid up to 30 kbar. For pressures exceeding this range, second-order pressure derivatives become significant and are not included in DensityX.

## Melt Composition

DensityX operates in a 10-component system: SiO{sub}`2`-TiO{sub}`2`-Al{sub}`2`O{sub}`3`-Fe{sub}`2`O{sub}`3`-FeO-MgO-CaO-Na{sub}`2`O−K{sub}`2`O-H{sub}`2`O.
**Figure 1** shows a Total-Alkali-Silica diagram showing the compositions of the calibration data.

```{figure} calibration_TAS.svg
:alt: TAS-diagram for calibration data
:align: center

**Figure 1.** Total-Alkali-Silica diagram for calibration data for studies used in DensityX. Kress1991, {cite:t}`kress1991`. Lange1997, {cite:t}`lange1997`. Ochs1999, {cite:t}`ochs1999`. Liu2006, {cite:t}`liu2006`. Guo2014, {cite:t}`guo2014`.
```

## Propagating Model Uncertainties

DensityX propagates the model parameter uncertainties in quadrature, where given, by assuming they are small, normally distributed and independent. The estimated uncertainty the calculated density is given as a {math}`1\sigma` value.

## References and Further Reading

```{bibliography} references.bib
:style: plain
