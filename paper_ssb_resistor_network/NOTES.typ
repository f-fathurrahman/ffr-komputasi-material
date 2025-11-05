
Title:
Using resistor network models to predict the transport properties of solid-state battery
composites

`https://doi.org/10.1038/s41467-025-56514-5`

Lukas Ketter, Niklas Greb, Tim Bernges and Wolfgang G. Zeier

---

Solid-state batteries utilize composite electrodes, consisting of
the electrochemically active material, the solid electrolyte (SE) and, if
necessary, additives such as binders or conductive carbons, e.g., vapor
grown carbon fibers (VGCF).

In these composite systems, a sufficiently
high ionic ($sigma_("ion")$) and electronic ($sigma_(e)$) conductivity are necessary, given
that both charge carriers need to access the active material during the
charging and discharging reactions.

Since a large mismatch in effective ionic and electronic conductivity
leads to inhomogeneous reaction
rates throughout the whole electrode thickness, balancing of the two
transport quantities is of paramount importance to enable high
electrode utilizations and avoid reaction fronts during battery operation.

---

On the one hand, charge carrier transport can be decisively
influenced by varying the volumetric ratios of the electrode
components given that the electrochemically active material and carbon
additives are solely responsible for the electronic conduction, while
the solid electrolyte exclusively conducts ions.

With that, altering the
SE content in composite anodes, such as Li4Ti5O12-Li7P2S8I-C10 or
Si-Li6PS5Cl (LPSCl)-C11, can facilitate changes in the effective ionic
conductivity over orders of magnitude.

Similarly strong correlations
between effective ionic and electronic conductivity and the volumetric
ratios of the composite components were found in cathode systems
such as NCM622-LPSCl12, LiMn2O4-Li3InCl6
13 and S-LPSCl-C14.

On the other hand, the microstructure of the composite electrode can
influence its effective conductivities.

Froboese et al. demonstrated the
influence of differently sized inclusions on the effective ionic
conductivity of composites.

Studies investigating the effect of changing active material
particle sizes in NCM622-LPSCl16 and Si-LPSCl-C17
electrodes on the battery performance highlight the importance of
carefully designing microstructures to improve electrode utilizations.

In a similar manner, optimizing the particle size of the SE in
Li Ni_0.83 Co_0.11 Mn_0.06 (NCM83)-LPSCl18 composites leads to better elec-
trode utilizations due to a more homogeneous ion current distribution
during battery operation.

This emphasizes that both the active material and
electrolyte particle size need to be considered in electrode
design.

---


While the influence of volume fractions and microstructure on
charge carrier transport is frequently investigated, the influence on the
thermal transport properties is scarcely studied.

Nonetheless, due to internal resistances ($R$) within the solid-state
battery system, Joule heating occurs during battery operation.

This is of particular importance during fast charging, as the power
of heating ($P$) is related to the square of the applied current ($I$)
via $P = R I^2$.

When the generation of
thermal energy is faster than its dissipation, heat accumulates in the
battery system leading to an increase of the battery temperature and
the development of temperature gradients.

This heating and temperature evolution within a cell can
significantly impact battery performance, accelerate degradation
kinetics and lead to at-times violent device failure.

Particularly considering battery safety, thermal
transport is known to play a critical role in conventional batteries.

Even though solid-state batteries are often promoted as being safer
than their conventional Li-ion counterparts, potential safety concerns
should not be overlooked and carefully studied.

Using thermodynamic models, Bates et al. showed, that under
short-circuit conditions larger temperature increases are possible
in SSB compared to batteries with liquid electrolytes.

Furthermore, severe reactions and fire generation
even under inert atmospheres were observed by Kim et al. when
they exposed charged cathode composites to elevated temperatures
above 150 °C.

This is further corroborated by exothermic reactions
that were observed between SSB-components in a differential scan-
ning calorimetry study by Johnson et al.

Especially regarding the
thermal management design, thermal transport considerations for SSB
pose an essential addition to the often-conducted ionic and electronic
conductivity studies.

Thermal transport investigations on various SE have shown that their
thermal conductivities are low and show a "glass-like" temperature-dependence.

Density dependent thermal transport
investigations on sulfidic $"Li"^(+)$ and $"Na"^(+)$ conducting
SE revealed further lowering of the thermal conductivity
with porosity, leading to thermal
conductivities in the range of liquid electrolytes.

The low thermal
conductivity of SE materials indicate that slow heat dissipation and
potential thermal runaway warrant careful consideration.

---

Generally, evaluating transport in composites as a function of
composition relies on a variety of measurements, making them a
time-consuming endeavor.

These measurements are tedious considering
that processing, additives, binders, and any microstructural changes
will influence the effective transport through the composite and that
transport needs to be re-evaluated every time.

Aside from experimental approaches, microstructure-based simulations
can offer valuable guidance when developing microstructural design concepts.

Considering the meso- and even atomistic scale, Heo et al.
investigated the impact of microstructure on ionic conductivity
in Li7 La3 Zr2 O12 computationally, while Haruyama et al.
calculated optimized structures and properties of interfaces
between oxidic cathode active material (CAM) and $beta$-Li3PS4
using density functional theory.

Microstructural modelling of entire composite cathodes can further
reveal percolation behavior of ion- and electron-conducting clusters,
while effective electronic and ionic conductivities can be calculated
using flux-based simulations on composite microstructure
representations.

Finsterbusch et al. constructed a virtual twin
based on scanning electron microscopy micrographs and
subsequently simulated discharge curves.

Extending this approach, Chen et al. utilized an electrochemical-thermal
coupled model to simulate heat evolution in addition to discharge
characteristics of composite cathodes.

Nevertheless, these approaches are often computationally
expensive, and a faster predictive approach is desirable to screen the
transport of new systems.
---

In this work, we propose a simple resistor network model
describing transport in solid-state battery electrode composites based
on simple voxel structures that is then benchmarked against
experimental transport in NCM83-LPSCl cathode composites.

A schematic
two-dimensional representation of a typical solid-state composite
cathode microstructure in comparison to a voxel microstructure
analogue is shown Fig. 1.

Despite this simplification, essential concepts
like tortuosity and domain size are reflected by the model.

Concluding the experimental characterization and microstructure-based transport
simulations of the ionic, electronic, and thermal transport properties
of the NCM83-LPSCl system, the proposed microstructural transport
model is tested against literature data to show its capability to predict
transport of other composite systems.

The results show that the model
can easily be modified and extended to other systems and research
questions.

Using the resistor network model provides experimentalists
a new tool to fast find the range of optimized solid-state battery
compositions - without the need for much computational resources.

---

= Results

== Resistor network modeling of solid-state battery composites

When a material is subjected to a temperature gradient ($nabla T$),
heat flow occurs.

According to Fourier's law (Eq. 1) the heat flux density ($q$) is
linked to the temperature gradient by the thermal conductivity ($kappa$)
of the material:
$
  q = -kappa nabla T
$
