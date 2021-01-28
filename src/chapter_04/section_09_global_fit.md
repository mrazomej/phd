Global Fit of All Parameters {#appendix_global_fit .unnumbered}
============================

In the main text, we used the repressor copy numbers $R$ and
repressor-DNA binding energies $\Delta\varepsilon_{RA}$ as reported by
@Garcia2011c. However, any error in these previous measurements of $R$
and $\Delta\varepsilon_{RA}$ will necessarily propagate into our own
fold-change predictions. In this section we take an alternative approach
to fitting the physical parameters of the system to that used in the
main text. First, rather than fitting only a single strain, we fit the
entire data set in along with microscopy data for the synthetic operator
Oid (see Appendix [\[AppendixOid\]](#AppendixOid){reference-type="ref"
reference="AppendixOid"}). In addition, we also simultaneously fit the
parameters $R$ and $\Delta\varepsilon_{RA}$ using the prior information
given by the previous measurements. By using the entire data set and
fitting all of the parameters, we obtain the best possible
characterization of the statistical mechanical parameters of the system
given our current state of knowledge. As a point of reference, we state
all of the parameters of the MWC model derived in the text in
Table [\[table_param\]](#table_param){reference-type="ref"
reference="table_param"}.

To fit all of the parameters simultaneously, we follow a similar
approach to the one detailed in the Methods section. Briefly, we perform
a Bayesian parameter estimation of the dissociation constants $K_A$ and
$K_I$, the six different repressor copy numbers $R$ corresponding to the
six *lacI* ribosomal binding sites used in our work, and the four
different binding energies $\Delta
\varepsilon_{RA}$ characterizing the four distinct operators used to
make the experimental strains. As in the main text, we fit the
logarithms $\tilde{k}_A =
-\log \frac{K_A}{1\,\text{M}}$ and
$\tilde{k}_I = -\log \frac{K_I}{1\,\text{M}}$ of the dissociation
constants which grants better numerical stability.

As in , we assume that deviations of the experimental fold-change from
the theoretical predictions are normally distributed with mean zero and
standard deviation $\sigma$. We begin by writing Bayes' theorem,
$$P(\tilde{k}_A, \tilde{k}_I, \mathbf{\textbf{\textit{R}}},
  \mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma \mid D) =
  \frac{P(D \mid \tilde{k}_A, \tilde{k}_I, \mathbf{\textbf{\textit{R}}},
  \mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma) P(\tilde{k}_A,
  \tilde{k}_I, \mathbf{\textbf{\textit{R}}},
  \mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma)}{P(D)},
  \label{eq_bayes_global}$$ where $\textbf{\textit{R}}$ is an array
containing the six different repressor copy numbers to be fit,
$\mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}$ is an
array containing the four binding energies to be fit, and $D$ is the
experimental fold-change data. The term $P(\tilde{k}_A, \tilde{k}_I,
\mathbf{\textbf{\textit{R}}},
\mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma \mid D)$
gives the probability distributions of all of the parameters given the
data. The term
$P(D \mid \tilde{k}_A, \tilde{k}_I, \mathbf{\textbf{\textit{R}}},
\mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma)$
represents the likelihood of having observed our experimental data given
some value for each parameter. $P(\tilde{k}_A, \tilde{k}_I,
\mathbf{\textbf{\textit{R}}},
\mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma)$
contains all the prior information on the values of these parameters.
Lastly, $P(D)$ serves as a normalization constant and hence can be
ignored.

Given $n$ independent measurements of the fold-change, the first term in
can be written as $$\begin{aligned}
    P(D \mid \tilde{k}_A, \tilde{k}_I, \mathbf{\textbf{\textit{R}}},
    \mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma) =
    \frac{1}{(2\pi\sigma^2)^{\frac{n}{2}}}\prod\limits_{i=1}^n \exp
    \left[-\frac{(\text{fc}^{(i)}_{\exp} - \text{fc}(\tilde{k}_A, \tilde{k}_I,
    R^{(i)}, \Delta\varepsilon_{RA}^{(i)}, c^{(i)}))^2}{2\sigma^2}\right],
    \label{eq_likelihood_global}\end{aligned}$$ where
$\text{fc}^{(i)}_{\text{exp}}$ is the $i^{\text{th}}$ experimental
fold-change and $\text{fc}(\cdot\cdot\cdot)$ is the theoretical
prediction. Note that the standard deviation $\sigma$ of this
distribution is not known and hence needs to be included as a parameter
to be fit.

The second term in represents the prior information of the parameter
values. We assume that all parameters are independent of each other, so
that $$P(\tilde{k}_A, \tilde{k}_I, \mathbf{\textbf{\textit{R}}},
\mathbf{\Delta\boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma) =
P(\tilde{k}_A ) \cdot P(\tilde{k}_I ) \cdot \prod_i P(R^{(i)}) \cdot \prod_j
P(\Delta\varepsilon_{RA}^{(j)}) \cdot P(\sigma),$$ where the superscript
$(i)$ indicates the repressor copy number of index $i$ and the
superscript $(j)$ denotes the binding energy of index $j$. As above, we
note that a prior must also be included for the unknown parameter
$\sigma$.

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
repressor copy numbers to be fit satisfies $$\label{eq_MCMC_R}
P(R^{(i)}) = \frac{1}{\sqrt{2\pi\sigma_{R_i}^2}} \exp \left(-\frac{(R^{(i)} -
    \bar{R}^{(i)})^2}{2 \sigma_{R_i}^2} \right),$$ where $\bar{R}^{(i)}$
is the mean repressor copy number and $\sigma_{R_i}$ is the variability
associated with this parameter as reported in @Garcia2011c. Note that we
use the given value of $\sigma_{R_i}$ from previous measurements rather
than leaving this as a free parameter.

Similarly, the binding energies $\Delta\varepsilon_{RA}^{(j)}$ are also
assumed to have a Gaussian informative prior of the same form. We write
it as
$$P(\Delta\varepsilon_{RA}^{(j)}) = \frac{1}{\sqrt{2\pi\sigma_{\varepsilon_j}^2}}
\exp \left(- \frac{(\Delta\varepsilon_{RA}^{(j)} -
    \Delta\bar{\varepsilon}_{RA}^{(j)})^2}{2 \sigma_{\varepsilon_j}^2} \right),$$
where $\Delta\bar{\varepsilon}_{RA}^{(j)}$ is the binding energy and
$\sigma_{\varepsilon_j}$ is the variability associated with that
parameter around the mean value as reported in @Garcia2011c .

The $\sigma_{R_i}$ and $\sigma_{\varepsilon_j}$ parameters will
constrain the range of values for $R^{(i)}$ and
$\Delta\varepsilon_{RA}^{(j)}$ found from the fitting. For example, if
for some $i$ the standard deviation $\sigma_{R_i}$ is very small, it
implies a strong confidence in the previously reported value.
Mathematically, the exponential in will ensure that the best-fit
$R^{(i)}$ lies within a few standard deviations of $\bar{R}^{(i)}$.
Since we are interested in exploring which values could give the best
fit, the errors are taken to be wide enough to allow the parameter
estimation to freely explore parameter space in the vicinity of the best
estimates. Putting all these terms together, we use Markov chain Monte
Carlo to sample the posterior distribution
$P(\tilde{k}_A, \tilde{k}_I, \mathbf{\textbf{\textit{R}}},
\mathbf{\Delta \boldsymbol{\varepsilon}_{\textbf{\textit{RA}}}}, \sigma \mid
D)$, enabling us to determine both the most likely value for each
physical parameter as well as its associated credible region (see the
[GitHub
repository](https://rpgroup-pboc.github.io/mwc_induction/code/notebooks/global_fits.html)
for the implementation).

shows the result of this global fit. When compared with we can see that
fitting for the binding energies and the repressor copy numbers improves
the agreement between the theory and the data.
Table [\[tab_global_fit\]](#tab_global_fit){reference-type="ref"
reference="tab_global_fit"} summarizes the values of the parameters as
obtained with this MCMC parameter inference. We note that even though we
allowed the repressor copy numbers and repressor-DNA binding energies to
vary, the resulting fit values were very close to the previously
reported values. The fit values of the repressor copy numbers were all
within one standard deviation of the previous reported values provided
in @Garcia2011c. And although some of the repressor-DNA binding energies
differed by a few standard deviations from the reported values, the
differences were always less than $1~k_BT$, which represents a small
change in the biological scales we are considering. The biggest
discrepancy between our fit values and the previous measurements arose
for the synthetic Oid operator, which we discuss in more detail in
Appendix [\[AppendixOid\]](#AppendixOid){reference-type="ref"
reference="AppendixOid"}.

shows the same key properties as in , but uses the parameters obtained
from this global fitting approach. We note that even by increasing the
number of degrees of freedom in our fit, the result does not change
substantially, due to in general, only minor improvements between the
theoretical curves and data. For the O3 operator data, again, agreement
between the predicted $[EC_{50}]$ and the effective Hill coefficient
remain poor due the theory being unable to capture the steepness of the
response curves.
