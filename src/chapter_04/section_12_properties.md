## Properties of Induction Titration Curves {#sec:ch4_sec12}

In this section, we expand on the phenotypic properties of the induction
response that were explored in the main text (see ). We begin by expanding on
our discussion of dynamic range and then show the analytic form of the
$[EC_{50}]$ for simple repression.

As stated in the main text, the dynamic range is defined as the difference
between the maximum and minimum system response, or equivalently, as the
difference between the saturation and leakiness of the system. Using , the
dynamic range is given by
$$
\text{dynamic range} = 
\left(
1+\frac{1}{1+e^{-\beta \Delta \varepsilon_{AI} }
 \left(\frac{K_A}{K_I}\right)^n }\frac{R}{N_{NS}}
 e^{-\beta \Delta\varepsilon_{RA}} \right)^{-1} - 
 \left(
1+\frac{1}{1+e^{-\beta \Delta \varepsilon_{AI} }}
\frac{R}{N_{NS}}e^{-\beta \Delta\varepsilon_{RA}} \right)^{-1}.
\label{eq:ch4_eq35}
$$
The dynamic range, along with saturation and leakiness were plotted with our
experimental data in - as a function of repressor copy number. shows how these
properties are expected to vary as a function of the repressor-operator binding
energy. Note that the resulting curves for all three properties have the same
shape as in -, since the dependence of the fold-change upon the repressor copy
number and repressor-operator binding energy are both contained in a single
multiplicative term, $R e^{-\beta \Delta\varepsilon_{RA}}$. Hence, increasing
$R$ on a logarithmic scale (as in -) is equivalent to decreasing
$\Delta\varepsilon_{RA}$ on a linear scale (as in ).

An interesting aspect of the dynamic range is that it exhibits a peak as a
function of either the repressor copy number (or equivalently of the
repressor-operator binding energy). Differentiating the dynamic range and
setting it equal to zero, we find that this peak occurs at
$$
\frac{R^*}{N_{NS}} = e^{-\beta (\Delta\varepsilon_{AI} -
\Delta\varepsilon_{RA})} \sqrt{e^{\Delta\varepsilon_{AI}} + 1} 
\sqrt{e^{\Delta\varepsilon_{AI}} + \left( \frac{K_A}{K_I} \right)^n}.
\label{eq:ch4_eq36}
$$
The magnitude of the peak is given by
$$
\text{max dynamic range} = 
\frac{\left( \sqrt{e^{\Delta\varepsilon_{AI}} + 1} - 
\sqrt{e^{\Delta\varepsilon_{AI}} + \left( \frac{K_A}{K_I} \right)^n} 
\right)^2}{\left( \frac{K_A}{K_I} \right)^n - 1},
\label{eq:ch4_eq37}
$$
which is independent of the repressor-operator binding energy
$\Delta\varepsilon_{RA}$ or $R$, and will only cause a shift in the location of
the peak but not its magnitude.

![**Dependence of leakiness, saturation, and dynamic range on the operator
binding energy and repressor copy number.** Increasing repressor copy number or
decreasing the repressor-operator binding energy suppresses gene expression and
decreases both the leakiness and saturation. The dynamic range retains its shape
but shifts right as the repressor copy number increases. The peak in the dynamic
range can be understood by considering the two extremes for $\Delta
\varepsilon_{RA}$: for small repressor-operator binding energies, the leakiness
is small but the saturation increases with $\Delta \varepsilon_{RA}$; for lare
repressor-operator binding energies the saturation is near unity and the
leakiness increases with $\Delta \varepsilon_{RA}$, thereby decreasing the
dynamic range. Repressor copy number does not affect the maximum dynamic range
(see ). Circles, diamonds, and squares represent $\Delta \varepsilon_{RA}$
values for the O1, O2, and O3 operators, respectively, demonstrating the
expected values of the properties using those
strains.](ch4_fig26){#fig:ch4_fig26 short-caption="Dependence of leakiness,
saturation, and dynamic range on the operator binding energy and repressor copy
number."}

We now consider the two remaining properties, the $[EC_{50}]$ and effective Hill
coefficient, which determine the horizontal properties of a system - that is,
they determine the range of inducer concentration in which the system's response
goes from its minimum to maximum values. The $[EC_{50}]$ denotes the inducer
concentration required to generate fold-change halfway between its minimum and
maximum value and was defined implicitly in . For the simple repression system,
the $[EC_{50}]$ is given by
$$
\frac{[EC_{50}]}{K_A} = 
\frac{\frac{K_A}{K_I} - 1}
{\frac{K_A}{K_I} - \left( \frac{\left( 1 + \frac{R}{N_{NS} }
e^{-\beta \Delta\varepsilon_{RA}} \right) + 
\left( \frac{K_A}{K_I} \right)^n 
\left( 2 e^{-\beta \Delta \varepsilon_{AI}} + 
\left( 1 + \frac{R}{N_{NS}} e^{-\beta \Delta\varepsilon_{RA}} \right) 
\right) }
{ 2 \left( 1 + \frac{R}{N_{NS}} e^{-\beta \Delta\varepsilon_{RA}} \right) +
 e^{-\beta \Delta \varepsilon_{AI}} + 
 \left( \frac{K_A}{K_I} \right)^n e^{-\beta \Delta \varepsilon_{AI}} } 
 \right)^{\frac{1}{n}}} - 1.
\label{eq:ch4_eq38}
$$
Using this expression, we can then find the effective Hill coefficient $h$,
which equals twice the log-log slope of the normalized fold-change evaluated at
$c = [EC_{50}]$ (see ). In - we show how these two properties vary with
repressor copy number, and in we demonstrate how they depend on the
repressor-operator binding energy. Both the $[EC_{50}]$ and $h$ vary
significantly with repressor copy number for sufficiently strong operator
binding energies. Interestingly, for weak operator binding energies on the order
of the O3 operator, it is predicted that the effective Hill coefficient should
not vary with repressor copy number. In addition, the maximum possible Hill
coefficient is roughly 1.75, which stresses the point that the effective Hill
coefficient should not be interpreted as the number of inducer binding sites,
which is exactly 2.

![**$\boldsymbol{[EC_{50}]}$ and effective Hill coefficient depend strongly on
repressor copy number and operator binding energy.** $[EC_{50}]$ values range
from very small and tightly clustered at weak operator binding energies (e.g.
O3) to relatively large and spread out for stronger operator binding energies
(O1 and O2). The effective Hill coefficient generally decreases with increasing
repressor copy number, indicating a flatter normalized response. The maximum
possible Hill coefficient is roughly 1.75 for all repressor-operator binding
energies. Circles, diamonds, and squares represent $\Delta \varepsilon_{RA}$
values for the O1, O2, and O3 operators,
respectively.](ch4_fig27){#fig:ch4_fig27
short-caption="**$\boldsymbol{[EC_{50}]}$ and effective Hill coefficient depend
strongly on repressor copy number and operator binding energy."}
