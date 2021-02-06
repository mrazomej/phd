## Empirical fits to noise predictions 

(Note: The Python code used for the calculations presented in this section can
be found in the [following
link](https://www.rpgroup.caltech.edu/chann_cap/src/theory/html/empirical_constants.html)
as an annotated Jupyter notebook)

In [Fig:ch3_fig03](C) in the main text we show that our minimal model has a
systematic deviation on the gene expression noise predictions compared to the
experimental data. This systematics will need to be addressed on an improved
version of the minimal model presented in this work. To guide the insights into
the origins of this systematic deviation in this appendix we will explore
empirical modifications of the model to improve the agreement between theory and
experiment.

### Multiplicative factor for the noise 

The first option we will explore is to modify our noise predictions by a
constant multiplicative factor. This means that we assume the relationship
between our minimal model predictions and the data for noise in gene expression
are of the from
\begin{equation}
\text{noise}_{\text{exp}} = \alpha \cdot \text{noise}_{\text{theory}},
\end{equation}
where $\alpha$ is a dimensionless constant to be fit from the data. The data,
especially in suggests that our predictions are within a factor of $\approx$ two
from the experimental data. To further check that intuition we performed a
weighted linear regression between the experimental and theoretical noise
measurements. The weight for each datum was taken to be proportional to the
bootstrap errors in the noise estimate, this to have poorly determined noises
weigh less during the regression. The result of this regression with no
intercept shows exactly that a factor of two systematically improves the
theoretical vs. experimental predictions. [@Fig:ch5_fig30] shows the improved
agreement when the theoretical predictions for the noise are multiplied by
$\approx 1.5$.

![**Multiplicative factor to improve theoretical vs. experimental comparison of
noise in gene expression.** Theoretical vs. experimental noise both in linear
(left) and log (right) scale. The dashed line shows the identity line of slope 1
and intercept zero. All data are colored by the corresponding value of the
experimental fold-change in gene expression as indicated by the color bar. The
$x$-axis was multiplied by a factor of $\approx 1.5$ as determined by a linear
regression from the data in . Each datum represents a single date measurement of
the corresponding strain and IPTG concentration with $\geq 300$ cells. The
points correspond to the median, and the error bars correspond to the 95%
confidence interval as determined by 10,000 bootstrap
samples.](ch5_fig30){#fig:ch5_fig30 short-caption="Multiplicative factor to
improve theoretical vs. experimental comparison of noise in gene expression"}

For completeness [@Fig:ch5_fig31] shows the noise in gene expression as a
function of the inducer concentration including this factor of $\approx 1.5$. It
is clear that overall a simple multiplicative factor improves the predictive
power of the model.

![**Protein noise of the regulated promoter with multiplicative factor.**
Comparison of the experimental noise for different operators ((A) O1,
$\Delta\varepsilon_r = -15.3 \; k_BT$, (B) O2, $\Delta\varepsilon_r = -13.9 \;
k_BT$, (C) O3, $\Delta\varepsilon_r = -9.7 \; k_BT$) with the theoretical
predictions for the the multi-promoter model. A linear regression revealed that
multiplying the theoretical noise prediction by a factor of $\approx 1.5$ would
improve agreement between theory and data. Points represent the experimental
noise as computed from single-cell fluorescence measurements of different *E.
coli* strains under 12 different inducer concentrations. Dotted line indicates
plot in linear rather than logarithmic scale. Each datum represents a single
date measurement of the corresponding strain and IPTG concentration with $\geq
300$ cells. The points correspond to the median, and the error bars correspond
to the 95% confidence interval as determined by 10,000 bootstrap samples.
White-filled dots are plot at a different scale for better
visualization.](ch5_fig31){#fig:ch5_fig31 short-caption="Protein noise of the
regulated promoter with multiplicative factor"}

### Additive factor for the noise 

As an alternative way to empirically improve the predictions of our model we
will now test the idea of an additive constant. What this means is that our
minimal model underestimates the noise in gene expression as
\begin{equation}
\text{noise}_{\text{exp}} = \beta + \text{noise}_{\text{theory}},
\end{equation}
where $\beta$ is an additive constant to be determined from the data. As with
the multiplicative constant we performed a regression to determine this
empirical additive constant comparing experimental and theoretical gene
expression noise values. We use the error in the 95% bootstrap confidence
interval as a weight for the linear regression. [@Fig:ch5_fig32] shows the
resulting theoretical vs. experimental noise where $\beta \approx 0.2$. We can
see a great improvement in the agreement between theory and experiment with this
additive constant.

![**Additive factor to improve theoretical vs. experimental comparison of noise
in gene expression.** Theoretical vs. experimental noise both in linear (left)
and log (right) scale. The dashed line shows the identity line of slope 1 and
intercept zero. All data are colored by the corresponding value of the
experimental fold-change in gene expression as indicated by the color bar. A
value of $\approx 0.2$ was added to all values in the $x$-axis as determined by
a linear regression from the data in . Each datum represents a single date
measurement of the corresponding strain and IPTG concentration with $\geq 300$
cells. The points correspond to the median, and the error bars correspond to the
95% confidence interval as determined by 10,000 bootstrap
samples.](ch5_fig32){#fig:ch5_fig32 short-caption="Additive factor to improve
theoretical vs. experimental comparison of noise in gene expression"}

For completeness [@Fig:ch5_fig33] shows the noise in gene expression as a
function of the inducer concentration including this additive factor of $\beta
\approx 0.2$. If anything, the additive factor seems to improve the agreement
between theory and data even more than the multiplicative factor.

![**Protein noise of the regulated promoter with additive factor.** Comparison
of the experimental noise for different operators ((A) O1, $\Delta\varepsilon_r
= -15.3 \; k_BT$, (B) O2, $\Delta\varepsilon_r = -13.9 \; k_BT$, (C) O3,
$\Delta\varepsilon_r = -9.7 \; k_BT$) with the theoretical predictions for the
the multi-promoter model. A linear regression revealed that an additive factor
of $\approx 0.2$ to the the theoretical noise prediction would improve agreement
between theory and data. Points represent the experimental noise as computed
from single-cell fluorescence measurements of different *E. coli* strains under
12 different inducer concentrations. Dotted line indicates plot in linear rather
than logarithmic scale. Each datum represents a single date measurement of the
corresponding strain and IPTG concentration with $\geq 300$ cells. The points
correspond to the median, and the error bars correspond to the 95% confidence
interval as determined by 10,000 bootstrap samples. White-filled dots are plot
at a different scale for better visualization.](ch5_fig33){#fig:ch5_fig33
short-caption="Protein noise of the regulated promoter with additive factor"}

### Correction factor for channel capacity with multiplicative factor

As seen in a constant multiplicative factor can reduce the discrepancy between
the model predictions and the data with respect to the noise (standard deviation
/ mean) in protein copy number. To find the equivalent correction would be for
the channel capacity requires gaining insights from the so-called small noise
approximation [@Tkacik2008a]. The small noise approximation assumes that the
input-output function can be modeled as a Gaussian distribution in which the
standard deviation is small. Using these assumptions one can derive a
closed-form for the channel capacity. Although our data and model predictions do
not satisfy the requirements for the small noise approximation, we can gain some
intuition for how the channel capacity would scale given a systematic deviation
in the cell-to-cell variability predictions compared with the data.

Using the small noise approximation one can derive the form of the input
distribution at channel capacity $P^*(c)$. To do this we use the fact that there
is a deterministic relationship between the input inducer concentration $c$ and
the mean output protein value $\left\langle p \right\rangle$, therefore we can
work with $P(\left\langle p \right\rangle)$ rather than $P(c)$ since the
deterministic relation allows us to write 
\begin{equation}
P(c) dc = P(\left\langle p \right\rangle) d\left\langle p \right\rangle.
\end{equation}
Optimizing over all possible distributions $P(\left\langle p \right\rangle)$
using calculus of variations results in a distribution of the form
\begin{equation}
P^*(\left\langle p \right\rangle) = 
\frac{1}{\mathcal{Z}} \frac{1}{\sigma_p(\left\langle p \right\rangle)},
\end{equation}
where $\sigma_p(\left\langle p \right\rangle)$ is the standard deviation of the
protein distribution as a function of the mean protein expression, and
$\mathcal{Z}$ is a normalization constant defined as
\begin{equation}
\mathcal{Z} \equiv 
\int_{\left\langle{p(c=0)}\right\rangle}
^{\left\langle{p(c\rightarrow \infty)}\right\rangle}
\frac{1}{\sigma_p(\left\langle p \right\rangle)} d\left\langle p \right\rangle.
\end{equation}
Under these assumptions the small noise approximation tells us that the channel
capacity is of the form [@Tkacik2008a]
\begin{equation}
I = \log_2 \left( \frac{\mathcal{Z}}{\sqrt{2 \pi e}} \right).
\end{equation}

From the theory-experiment comparison in we know that the standard deviation
predicted by our model is systematically off by a factor of two compared to the
experimental data, i.e.
\begin{equation}
\sigma_p^{\exp} = 2 \sigma_p^{\text{theory}}.
\end{equation}
This then implies that the normalization constant $\mathcal{Z}$ between theory
and experiment must follow a relationship of the form
\begin{equation}
\mathcal{Z}^{\exp} = \frac{1}{2} \mathcal{Z}^{\text{theory}}.
\end{equation}
With this relationship the small noise approximation would predict that the
difference between the experimental and theoretical channel capacity should be
of the form
\begin{equation}
I^{\exp} = \log_2 \left( \frac{\mathcal{Z}^{\exp}}{\sqrt{2 \pi e}} \right)
= \log_2 \left( \frac{\mathcal{Z}^{\text{theory}}}{\sqrt{2 \pi e}} \right)
- \log_2(2).
\end{equation}
Therefore under the small noise approximation we would expect our predictions
for the channel capacity to be off by a constant of 1 bit ($\log_2(2)$) of
information. Again, the conditions for the small noise approximation do not
apply to our data given the intrinsic level of cell-to-cell variability in the
system, nevertheless what this analysis tells is is that we expect that an
additive constant should be able to explain the discrepancy between our model
predictions and the experimental channel capacity. To test this hypothesis we
performed a "linear regression" between the model predictions and the
experimental channel capacity with a fixed slope of 1. The intercept of this
regression, -0.56 bits, indicates the systematic deviation we expect should
explain the difference between our model and the data. [@Fig:ch5_fig34] shows
the comparison between the original predictions shown in (A) and the resulting
predictions with this shift. Other than the data with zero channel capacity,
this shift is able to correct the systematic deviation for all data. We
therefore conclude that our model ends up underestimating the experimentally
determined channel capacity by a constant amount of 0.43 bits.

![**Additive correction factor for channel capacity.** Solid lines represent the
theoretical predictions of the channel capacity shown in (A). The dashed lines
show the resulting predictions with a constant shift of -0.43 bits. Points
represent single biological replicas of the inferred channel
capacity.](ch5_fig34){#fig:ch5_fig34 short-caption="Additive correction factor
for channel capacity"}