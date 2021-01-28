## Global Fit of All Parameters

In the main text, we used the repressor copy numbers $R$ and repressor-DNA
binding energies $\Delta\varepsilon_{RA}$ as reported by @Garcia2011c. However,
any error in these previous measurements of $R$ and $\Delta\varepsilon_{RA}$
will necessarily propagate into our own fold-change predictions. In this section
we take an alternative approach to fitting the physical parameters of the system
to that used in the main text. First, rather than fitting only a single strain,
we fit the entire data set in along with microscopy data for the synthetic
operator Oid (see Appendix XXX). In addition, we also simultaneously fit the
parameters $R$ and $\Delta\varepsilon_{RA}$ using the prior information given by
the previous measurements. By using the entire data set and fitting all of the
parameters, we obtain the best possible characterization of the statistical
mechanical parameters of the system given our current state of knowledge. As a
point of reference, we state all of the parameters of the MWC model derived in
the text in Table XXX.

To fit all of the parameters simultaneously, we follow a similar approach to the
one detailed in the Methods section. Briefly, we perform a Bayesian parameter
estimation of the dissociation constants $K_A$ and $K_I$, the six different
repressor copy numbers $R$ corresponding to the six *lacI* ribosomal binding
sites used in our work, and the four different binding energies $\Delta
\varepsilon_{RA}$ characterizing the four distinct operators used to make the
experimental strains. As in the main text, we fit the logarithms $\tilde{k}_A =
-\log \frac{K_A}{1\,\text{M}}$ and $\tilde{k}_I = -\log \frac{K_I}{1\,\text{M}}$
of the dissociation constants which grants better numerical stability.

As in , we assume that deviations of the experimental fold-change from the
theoretical predictions are normally distributed with mean zero and standard
deviation $\sigma$. We begin by writing Bayes' theorem,
$$
P(\tilde{k}_A, \tilde{k}_I, \mathbf{\textbf{\textit{R}}},
\mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma \mid D) =
\frac{P(D \mid \tilde{k}_A, \tilde{k}_I, \mathbf{\textbf{\textit{R}}},
\mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma) P(\tilde{k}_A,
\tilde{k}_I, \mathbf{\textbf{\textit{R}}},
\mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma)}{P(D)},
$${#eq:ch4_eq30}
where $\textbf{\textit{R}}$ is an array containing the six different repressor
copy numbers to be fit,
$\mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}$ is an array
containing the four binding energies to be fit, and $D$ is the experimental
fold-change data. The term $P(\tilde{k}_A, \tilde{k}_I,
\mathbf{\textbf{\textit{R}}},
\mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma \mid D)$
gives the probability distributions of all of the parameters given the data. The
term $P(D \mid \tilde{k}_A, \tilde{k}_I, \mathbf{\textbf{\textit{R}}},
\mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma)$
represents the likelihood of having observed our experimental data given some
value for each parameter. $P(\tilde{k}_A, \tilde{k}_I,
\mathbf{\textbf{\textit{R}}},
\mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma)$
contains all the prior information on the values of these parameters. Lastly,
$P(D)$ serves as a normalization constant and hence can be
ignored.

Given $n$ independent measurements of the fold-change, the first term in
can be written as 
$$\begin{split}
P(D \mid \tilde{k}_A, \tilde{k}_I, \mathbf{\textbf{\textit{R}}},
\mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma) =
\frac{1}{(2\pi\sigma^2)^{\frac{n}{2}}}\prod\limits_{i=1}^n \exp
\left[-\frac{(\text{fc}^{(i)}_{\exp} - \text{fc}(\tilde{k}_A, \tilde{k}_I,
R^{(i)}, \Delta\varepsilon_{RA}^{(i)}, c^{(i)}))^2}{2\sigma^2}\right],
\end{split}
$${#eq:ch4_eq31}
where $\text{fc}^{(i)}_{\text{exp}}$ is the $i^{\text{th}}$ experimental
fold-change and $\text{fc}(\cdot\cdot\cdot)$ is the theoretical prediction. Note
that the standard deviation $\sigma$ of this distribution is not known and hence
needs to be included as a parameter to be fit.

The second term in represents the prior information of the parameter
values. We assume that all parameters are independent of each other, so
that 
$$
P(\tilde{k}_A, \tilde{k}_I, \mathbf{\textbf{\textit{R}}},
\mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma) =
P(\tilde{k}_A ) \cdot P(\tilde{k}_I ) \cdot \prod_i P(R^{(i)}) \cdot \prod_j
P(\Delta\varepsilon_{RA}^{(j)}) \cdot P(\sigma),
$${#eq:ch4_eq32}
where the superscript $(i)$ indicates the repressor copy number of index $i$ and
the superscript $(j)$ denotes the binding energy of index $j$. As above, we note
that a prior must also be included for the unknown parameter $\sigma$.

Because we knew nothing about the values of $\tilde{k}_A$,
$\tilde{k}_I$, and $\sigma$ before performing the experiment, we assign
maximally uninformative priors to each of these parameters. More
specifically, we assign uniform priors to $\tilde{k}_A$ and
$\tilde{k}_I$ and a Jeffreys prior to $\sigma$, indicating that $K_A$,
$K_I$, and $\sigma$ are scale parameters [@Sivia2006]. We do, however,
have prior information for the repressor copy numbers and the
repressor-DNA binding energies from @Garcia2011c. This prior knowledge is
included within our model using an informative prior for these two
parameters, which we assume to be Gaussian. Hence each of the $R^{(i)}$
repressor copy numbers to be fit satisfies
$$
P(R^{(i)}) = \frac{1}{\sqrt{2\pi\sigma_{R_i}^2}} 
\exp \left(-\frac{(R^{(i)} - \bar{R}^{(i)})^2}{2 \sigma_{R_i}^2} \right),
$${#eq:ch4_eq33}
where $\bar{R}^{(i)}$ is the mean repressor copy number and $\sigma_{R_i}$ is
the variability associated with this parameter as reported in @Garcia2011c. Note
that we use the given value of $\sigma_{R_i}$ from previous measurements rather
than leaving this as a free parameter.

Similarly, the binding energies $\Delta\varepsilon_{RA}^{(j)}$ are also assumed
to have a Gaussian informative prior of the same form. We write it as
$$
P(\Delta\varepsilon_{RA}^{(j)}) = 
\frac{1}{\sqrt{2\pi\sigma_{\varepsilon_j}^2}}
\exp \left(- \frac{(\Delta\varepsilon_{RA}^{(j)} -
\Delta\bar{\varepsilon}_{RA}^{(j)})^2}{2 \sigma_{\varepsilon_j}^2} \right),
$${#eq:ch4_eq34}
where $\Delta\bar{\varepsilon}_{RA}^{(j)}$ is the binding energy and
$\sigma_{\varepsilon_j}$ is the variability associated with that
parameter around the mean value as reported in @Garcia2011c .

The $\sigma_{R_i}$ and $\sigma_{\varepsilon_j}$ parameters will constrain the
range of values for $R^{(i)}$ and $\Delta\varepsilon_{RA}^{(j)}$ found from the
fitting. For example, if for some $i$ the standard deviation $\sigma_{R_i}$ is
very small, it implies a strong confidence in the previously reported value.
Mathematically, the exponential in will ensure that the best-fit $R^{(i)}$ lies
within a few standard deviations of $\bar{R}^{(i)}$. Since we are interested in
exploring which values could give the best fit, the errors are taken to be wide
enough to allow the parameter estimation to freely explore parameter space in
the vicinity of the best estimates. Putting all these terms together, we use
Markov chain Monte Carlo to sample the posterior distribution $P(\tilde{k}_A,
\tilde{k}_I, \mathbf{\textbf{\textit{R}}}, \mathbf{\Delta
\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma \mid D)$, enabling us
to determine both the most likely value for each physical parameter as well as
its associated credible region (see the [GitHub
repository](https://rpgroup-pboc.github.io/mwc_induction/code/notebooks/global_fits.html)
for the implementation).

shows the result of this global fit. When compared with we can see that fitting
for the binding energies and the repressor copy numbers improves the agreement
between the theory and the data. Table XXX summarizes the values of the
parameters as obtained with this MCMC parameter inference. We note that even
though we allowed the repressor copy numbers and repressor-DNA binding energies
to vary, the resulting fit values were very close to the previously reported
values. The fit values of the repressor copy numbers were all within one
standard deviation of the previous reported values provided in @Garcia2011c. And
although some of the repressor-DNA binding energies differed by a few standard
deviations from the reported values, the differences were always less than
$1~k_BT$, which represents a small change in the biological scales we are
considering. The biggest discrepancy between our fit values and the previous
measurements arose for the synthetic Oid operator, which we discuss in more
detail in Appendix XXX.

XXX shows the same key properties as in XXX, but uses the parameters obtained
from this global fitting approach. We note that even by increasing the number of
degrees of freedom in our fit, the result does not change substantially, due to
in general, only minor improvements between the theoretical curves and data. For
the O3 operator data, again, agreement between the predicted $[EC_{50}]$ and the
effective Hill coefficient remain poor due the theory being unable to capture
the steepness of the response curves.

| \textbf{Parameter}                               | \multicolumn{1}{c}{\textbf{Description}}                                                                     |
| ------------------------------------------------ | ------------------------------------------------------------------------------------------------------------ |
| $c$                                              | Concentration of the inducer                                                                                 |
| $K_A, K_I$                                       | Dissociation constant between an inducer and the repressor in the active/inactive state                      |
| $\Delta \varepsilon_{AI}$                        | The difference between the free energy of repressor in the inactive and active states                        |
| $\Delta\varepsilon_{P}$                          | Binding energy between the RNAP and its specific binding site                                                |
| $\Delta\varepsilon_{RA}, \Delta\varepsilon_{RI}$ | Binding energy between the operator and the active/inactive repressor                                        |
| $n$                                              | Number of inducer binding sites per repressor                                                                |
| $P$                                              | Number of RNAP                                                                                               |
| $R_A, R_I, R$                                    | Number of active/inactive/total repressors                                                                   |
| $p_A = \frac{R_A}{R}$                            | Probability that a repressor will be in the active state                                                     |
| $p_{\text{bound}}$                               | Probability that an RNAP is bound to the promoter of interest, assumed to be proportional to gene expression |
| $\text{fold-change}$                             | Ratio of gene expression in the presence of repressor to that in the absence of repressor                    |
| $F$                                              | Free energy of the system                                                                                    |
| $N_{NS}$                                         | The number of non-specific binding sites for the repressor in the genome                                     |
| $\beta = \frac{1}{k_B T}$                        | The inverse product of the Boltzmann constant $k_B$ and the temperature $T$ of the system                    |
Table: **Key model parameters for induction of an allosteric repressor.** {#tbl:ch4_tbl02}

![**Global fit of dissociation constants, repressor copy numbers and binding
energies.** Theoretical predictions resulting from simultaneously fitting the
dissociation constants $K_A$ and $K_I$, the six repressor copy numbers $R$, and
the four repressor-DNA binding energies $\Delta\varepsilon_{RA}$ using the
entire data set from XXX as well as the microscopy data for the Oid operator.
Error bars of experimental data show the standard error of the mean (eight or
more replicates) and shaded regions denote the 95\% credible region. Where error
bars are not visible, they are smaller than the point itself. For the Oid
operator, all of the data points are shown since a smaller number of replicates
were taken. The shaded regions are significantly smaller than in XXX because
this fit was based on all data points, and hence the fit parameters are much
more tightly constrained. The dashed lines at 0 IPTG indicates a linear scale,
whereas solid lines represent a log scale.](ch4_fig19){#fig:ch4_fig19
short-caption="Global fit of dissociation constants, repressor copy numbers and
binding energies"}

![**Key properties of induction profiles as predicted with a global fit using
all available data.** Data for the (A) leakiness, (B) saturation, and (C)
dynamic range are obtained from fold-change measurements in XXX in the absence
and presence of IPTG. All prediction curves were generated using the parameters
listed in \ref{tab_global_fit}. Both the (D) $[EC_{50}]$ and (E) effective Hill
coefficient are inferred by individually fitting all parameters--$K_A,\, K_I,\,
R,\, \Delta\varepsilon_{RA}$--to each operator-repressor pairing in XXX(A)-(C)
separately to [@Eq:fold_change_full] in order to smoothly interpolate between
the data points. Note that where error bars are not visible, this indicates that
the error bars are smaller than the point itself.](ch4_fig20){#fig:ch4_fig20
short-caption="Key properties of induction profiles as predicted with a global
fit using all available data"}

|                              | \textbf{Reported Values} [@Garcia2011c] | \textbf{Global Fit}                 |
| ---------------------------- | --------------------------------------- | ----------------------------------- |
| $\tilde{k}_A$                | $-$                                     | $-5.33^{+0.06}_{-0.05}$             |
| $\tilde{k}_I$                | $-$                                     | $0.31^{+0.05}_{-0.06}$              |
| $K_A$                        | $-$                                     | $205^{+11}_{-12}\,\mu\text{M}$      |
| $K_I$                        | $-$                                     | $0.73^{+0.04}_{-0.04}\,\mu\text{M}$ |
| $R_{22}$                     | $22 \pm 4$                              | $20^{+1}_{-1}$                      |
| $R_{60}$                     | $60 \pm 20$                             | $74^{+4}_{-3}$                      |
| $R_{124}$                    | $124 \pm 30$                            | $130^{+6}_{-6}$                     |
| $R_{260}$                    | $260 \pm 40$                            | $257^{+9}_{-11}$                    |
| $R_{1220}$                   | $1220 \pm 160$                          | $1191^{+32}_{-55}$                  |
| $R_{1740}$                   | $1740 \pm 340$                          | $1599^{+75}_{-87}$                  |
| O1 $\Delta\varepsilon_{RA}$  | $-15.3 \pm 0.2~k_BT$                    | $-15.2^{+0.1}_{-0.1}~k_BT$          |
| O2 $\Delta\varepsilon_{RA}$  | $-13.9 \pm 0.2~k_BT$                    | $-13.6^{+0.1}_{-0.1}~k_BT$          |
| O3 $\Delta\varepsilon_{RA}$  | $-9.7 \pm 0.1~k_BT$                     | $-9.4^{+0.1}_{-0.1}~k_BT$           |
| Oid $\Delta\varepsilon_{RA}$ | $-17.0 \pm 0.2~k_BT$                    | $-17.7^{+0.2}_{-0.1}~k_BT$          |
Table: **Global fit of all parameter values using the entire data set in XXX.**
In addition to fitting the repressor inducer dissociation constants $K_A$ and
$K_I$ as was done in the text, we also fit the repressor DNA binding energy
$\Delta\varepsilon_{RA}$ as well as the repressor copy numbers $R$ for each
strain. The middle columns show the previously reported values for all
$\Delta\varepsilon_{RA}$ and $R$ values, with $\pm$ representing the standard
deviation of three replicates. The right column shows the global fits from this
work, with the subscript and superscript notation denoting the 95\% credible
region. Note that there is overlap between all of the repressor copy numbers and
that the net difference in the repressor-DNA binding energies is less than
$1~k_B T$. The logarithms $\tilde{k}_A = -\log \frac{K_A}{1\,\text{M}}$ and
$\tilde{k}_I = -\log \frac{K_I}{1\,\text{M}}$ of the dissociation constants were
fit for numerical stability. {#tbl:ch4_tbl03}
