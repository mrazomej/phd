Applicability of Theory to the Oid Operator Sequence {#AppendixOid}
====================================================

In addition to the native operator sequences (O1, O2, and O3) considered
in the main text, we were also interested in testing our model
predictions against the synthetic Oid operator. In contrast to the other
operators, Oid is one base pair shorter in length ($20\,\text{bp}$), is
fully symmetric, and is known to provide stronger repression than the
native operator sequences considered so far. While the theory should be
similarly applicable, measuring the lower fold-changes associated with
this YFP construct was expected to be near the sensitivity limit for our
flow cytometer, due to the especially strong binding energy of Oid
($\Delta
\varepsilon_{RA}=-17.0 ~k_BT$) [@Garcia2011b]. Accordingly, fluorescence
data for Oid were obtained using microscopy, which is more sensitive
than flow cytometry. Appendix
[\[AppendixMicroscopy\]](#AppendixMicroscopy){reference-type="ref"
reference="AppendixMicroscopy"} gives a detailed explanation of how
microscopy measurements were used to obtain induction curves.

We follow the approach of the main text and make fold-change predictions
based on the parameter estimates from our strain with $R=260$ and an O2
operator. These predictions are shown in , where we also plot data taken
in triplicate for strains containing $R= 22$, 60, and 124, obtained by
single-cell microscopy. We find that the data are systematically below
the theoretical predictions. We also considered our global fitting
approach (see Appendix
[\[appendix_global_fit\]](#appendix_global_fit){reference-type="ref"
reference="appendix_global_fit"}) to see whether we might find better
agreement with the observed data. Interestingly, we findthat the
majority of the parameters remain largely unchanged, but our estimate
for the Oid binding energy $\Delta \varepsilon_{RA}$ is shifted to
$-17.7~k_BT$ instead of the value $-17.0~k_BT$ found by @Garcia2011c. In
we again plot the Oid fold-change data but with theoretical predictions
using the new estimate for the Oid binding energy from our global fit
and find substantially better agreement.

![**Predictions of fold-change for strains with an Oid binding sequence
versus experimental measurements with different repressor copy
numbers.** Experimental data is plotted against the parameter-free
predictions that are based on our fit to the O2 strain with $R=260$.
Here we use the previously measured binding energy
$\Delta\varepsilon_{RA}=-17.0~k_BT$ [@Garcia2011c]. The same experimental
data is plotted against the best-fit parameters using the complete O1,
O2, O3, and Oid data sets to infer $K_A$, $K_I$, repressor copy numbers,
and the binding energies of all operators (see Appendix
[\[appendix_global_fit\]](#appendix_global_fit){reference-type="ref"
reference="appendix_global_fit"}). Here the major difference in the
inferred parameters is a shift in the binding energy for Oid from
$\Delta\varepsilon_{RA}=-17.0~k_BT$ to
$\Delta\varepsilon_{RA}=-17.7~k_BT$, which now shows agreement between
the theoretical predictions and experimental data. Shaded regions from
the theoretical curves denote the 95% credible region. These are
narrower in Panel because the inference of parameters was performed with
much more data, and hence the best-fit values are more tightly
constrained. Individual data points are shown due to the small number of
replicates. The dashed lines at 0 IPTG indicate a linear scale, whereas
solid lines represent a log
scale.](SI_figs/figS21.pdf){#fig_Oid_theory_compare}

shows the cumulative data from @Garcia2011c and @Brewster2014, as well as
our data with $c=0 \, \mu \text{M}$, which all measured fold-change for
the same simple repression architecture utilizing different reporters
and measurement techniques. We find that the binding energies from the
global fit, including $\Delta \varepsilon_{RA}=-17.7~k_BT$, compare
reasonably well with all previous measurements.

![**Comparison of fold-change predictions based on binding energies from
Garcia and Phillips and those inferred from this work.** Fold-change
curves for the different repressor-DNA binding energies
$\Delta\varepsilon_{RA}$ are plotted as a function of repressor copy
number when IPTG concentration $c=0$. Solid curves use the binding
energies determined from @Garcia2011c, while the dashed curves use the
inferred binding energies we obtained when performing a global fit of
$K_A$, $K_I$, repressor copy numbers, and the binding energies using all
available data from our work. Fold-change measurements from our
experiments (outlined circles) @Garcia2011c (solid circles), and
@Brewster2014 (diamonds) show that the small shifts in binding energy
that we infer are still in agreement with prior data. Note that only a
single flow cytometry data point is shown for Oid from this study, since
the $R=60$ and $R=124$ curves from had extremely low fold-change in the
absence of inducer ($c=0$) so as to be indistinguishable from
autofluorescence, and in fact their fold-change values in this limit
were negative and hence do not appear on this
plot.](SI_figs/figS22.pdf){#fig_titration_summary_new_Oid_energy}
