## Comparison of Parameter Estimation and Fold-Change Predictions across Strains {#sec:ch4_sec11}

The inferred parameter values for $K_A$ and $K_I$ in the main text were
determined by fitting to induction fold-change measurements from a single strain
($R=260$, $\Delta\varepsilon_{RA} = -13.9~k_BT$, $n=2$, and
$\Delta\varepsilon_{AI}=4.5~k_BT$). After determining these parameters, we were
able to predict the fold-change of the remaining strains without any additional
fitting. However, the theory should be independent of the specific strain used
to estimate $K_A$ and $K_I$; using any alternative strain to fit $K_A$ and $K_I$
should yield similar predictions. For the sake of completeness, here we discuss
the values for $K_A$ and $K_I$ that are obtained by fitting to each of the
induction data sets individually. These fit parameters are shown in
[@Fig:ch2_fig04](D), where we find close agreement between strains, but with
some deviation and poorer inferences observed with the O3 operator strains.
Overall, we find that regardless of which strain is chosen to determine the
unknown parameters, the predictions laid out by the theory closely match the
experimental measurements. Here we present a comparison of the strain-specific
predictions and measured fold-change data for each of the three operators
considered.

We follow the approach taken in the main text and use Eq.
$\ref{eq:fold_change_full}$ to infer values for $K_A$ and $K_I$ by fitting to
each combination of binding energy $\Delta \varepsilon_{RA}$ and repressor copy
number $R$. We then use these fitted parameters to predict the induction curves
of all other strains. In [@Fig:ch4_fig23] we plot these fold-change predictions
along with experimental data for each of our strains that contain an O1
operator. To make sense of this plot, consider the first row as an example. In
the first row, $K_A$ and $K_I$ were estimated using data from the strain
containing $R=22$ and an O1 operator (top leftmost plot, shaded in gray). The
remaining plots in this row show the predicted fold-change using these values
for $K_A$ and $K_I$. In each row, we then infer $K_A$ and $K_I$ using data from
a strain containing a different repressor copy number ($R=60$ in the second row,
$R=124$ in the third row, and so on). In [@Fig:ch4_fig24] and [@Fig:ch4_fig25],
we similarly apply this inference to our strains with O2 and O3 operators,
respectively. We note that the overwhelming majority of predictions closely
match the experimental data. The notable exception is that using the $R=22$
strain provides poor predictions for the strains with large copy numbers
(especially $R=1220$ and $R=1740$), though it should be noted that predictions
made from the $R=22$ strain have considerably broader credible regions. This
loss in predictive power is due to the poorer estimates of $K_A$ and $K_I$ for
the $R=22$ strain as shown in [@Fig:ch2_fig04](D).

![**O1 strain fold-change predictions based on strain-specific parameter
estimation of $\boldsymbol{K_A}$ and $\boldsymbol{K_I}$.** Fold-change in
expression is plotted as a function of IPTG concentration for all strains
containing an O1 operator. The solid points correspond to the mean experimental
value. The solid lines correspond to Eq. $\ref{eq:fold_change_full}$ using the
parameter estimates of $K_A$ and $K_I$. Each row uses a single set of parameter
values based on the strain noted on the left axis. The shaded plots along the
diagonal are those where the parameter estimates are plotted along with the data
used to infer them. Values for repressor copy number and operator binding energy
are from [@Garcia2011c]. The shaded region on the curve represents the
uncertainty from our parameter estimates and reflects the 95% highest
probability density region of the parameter
predictions.](ch4_fig23){#fig:ch4_fig23 short-caption="O1 strain fold-change
predictions based on strain-specific parameter estimation of $\boldsymbol{K_A}$
and $\boldsymbol{K_I}$"}

![**O2 strain fold-change predictions based on strain-specific parameter
estimation of $\boldsymbol{K_A}$ and $\boldsymbol{K_I}$.** Fold-change in
expression is plotted as a function of IPTG concentration for all strains
containing an O2 operator. The plots and data shown are analogous to
[@Fig:ch4_fig23], but for the O2 operator.](ch4_fig24){#fig:ch4_fig24
short-caption="O2 strain fold-change predictions based on strain-specific
parameter estimation of $\boldsymbol{K_A}$ and $\boldsymbol{K_I}$"}

![**O3 strain fold-change predictions based on strain-specific parameter
estimation of $\boldsymbol{K_A}$ and $\boldsymbol{K_I}$.** Fold-change in
expression is plotted as a function of IPTG concentration for all strains
containing an O3 operator. The plots and data shown are analogous to
[@Fig:ch4_fig23], but for the O3 operator. We note that when using the $R=22$ O3
strain to predict $K_A$ and $K_I$, the large uncertainty in the estimates of
these parameters (see [@Fig:ch2_fig04](D)) leads to correspondingly wider
credible regions.](ch4_fig25){#fig:ch4_fig25 short-caption="O3 strain
fold-change predictions based on strain-specific parameter estimation of
$\boldsymbol{K_A}$ and $\boldsymbol{K_I}$"}
