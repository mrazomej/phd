## Three-state promoter model for simple repression {#supp_model}

In order to tackle the question of how much information the simple repression
motif can process we require the joint probability distribution of mRNA and
protein $P(m, p; t)$. To obtain this distribution we use the chemical master
equation formalism as described in . Specifically, we assume a three-state model
where the promoter can be found 1) in a transcriptionally active state ($A$
state), 2) in a transcriptionally inactive state without the repressor bound
($I$ state) and 3) with the repressor bound ($R$ state). (See (A)). These three
states generate a system of coupled differential equations for each of the three
state distributions $P_A(m, p)$, $P_I(m, p)$ and $P_R(m, p)$. Given the rates
shown in (A) let us define the system of ODEs. For the transcriptionally active
state we have 
\begin{equation}
\begin{split}
    \frac{d P_A(m, p)}{dt} &=
    - \overbrace{k^{(p)}_{\text{off}} P_A(m, p)}^{A \rightarrow I} % A -> I
    + \overbrace{k^{(p)}_{\text{on}} P_I(m, p)}^{I \rightarrow A}\\ % I -> A
    &+ \overbrace{r_m P_A(m-1, p)}^{m-1 \rightarrow m} % m-1 -> m
    - \overbrace{r_m P_A(m, p)}^{m \rightarrow m+1}% m -> m+1
    + \overbrace{\gamma _m (m + 1) P_A(m+1 , p)}^{m+1 \rightarrow m} % m+1 -> m
    - \overbrace{\gamma _m m P_A(m , p)}^{m \rightarrow m-1}\\ % m -> m-1
    &+ \overbrace{r_p m P_A(m, p - 1)}^{p-1 \rightarrow p} % p-1 -> p
    - \overbrace{r_p m P_A(m, p)}^{p \rightarrow p+1} % p -> p+1
    + \overbrace{\gamma _p (p + 1) P_A(m, p + 1)}^{p + 1 \rightarrow p} % p+1 -> p
    - \overbrace{\gamma _p p P_A(m, p)}^{p \rightarrow p-1}. % p -> p-1
\end{split}
\end{equation}
For the inactive promoter state $I$ we have
\begin{equation}
\begin{split}
    \frac{d P_I(m, p)}{dt} &=
    \overbrace{k^{(p)}_{\text{off}} P_A(m, p)}^{A \rightarrow I} % A -> I
    - \overbrace{k^{(p)}_{\text{on}} P_I(m, p)}^{I \rightarrow A} % I -> A
    + \overbrace{k^{(r)}_{\text{off}} P_R(m, p)}^{R \rightarrow I} % R -> I
    - \overbrace{k^{(r)}_{\text{on}} P_I(m, p)}^{I \rightarrow R}\\ % I -> R
    &+ \overbrace{\gamma _m (m + 1) P_I(m+1 , p)}^{m+1 \rightarrow m} % m+1 -> m
    - \overbrace{\gamma _m m P_I(m , p)}^{m \rightarrow m-1}\\ % m -> m-1
    &+ \overbrace{r_p m P_I(m, p - 1)}^{p-1 \rightarrow p} % p-1 -> p
    - \overbrace{r_p m P_I(m, p)}^{p \rightarrow p+1} % p -> p+1
    + \overbrace{\gamma _p (p + 1) P_I(m, p + 1)}^{p + 1 \rightarrow p} % p+1 -> p
    - \overbrace{\gamma _p p P_I(m, p)}^{p \rightarrow p-1}. % p -> p-1
  \end{split}
\end{equation}
And finally for the repressor bound state $R$ we have
\begin{equation}
\begin{split}
    \frac{d P_R(m, p)}{dt} &=
    - \overbrace{k^{(r)}_{\text{off}} P_R(m, p)}^{R \rightarrow I} % R -> I
    + \overbrace{k^{(r)}_{\text{on}} P_I(m, p)}^{I \rightarrow R}\\ % I -> R
    &+ \overbrace{\gamma _m (m + 1) P_R(m+1 , p)}^{m+1 \rightarrow m} % m+1 -> m
    - \overbrace{\gamma _m m P_R(m , p)}^{m \rightarrow m-1}\\ % m -> m-1
    &+ \overbrace{r_p m P_R(m, p - 1)}^{p-1 \rightarrow p} % p-1 -> p
    - \overbrace{r_p m P_R(m, p)}^{p \rightarrow p+1} % p -> p+1
    + \overbrace{\gamma _p (p + 1) P_R(m, p + 1)}^{p + 1 \rightarrow p} % p+1 -> p
    - \overbrace{\gamma _p p P_R(m, p)}^{p \rightarrow p-1}. % p -> p-1
\end{split}
\end{equation}

For an unregulated promoter, i.e. a promoter in a cell that has no repressors
present, and therefore constitutively expresses the gene, we use a two-state
model in which the state $R$ is not allowed. All the terms in the system of ODEs
containing $k^{(r)}_{\text{on}}$ or $k^{(r)}_{\text{off}}$ are then set to zero.

As detailed in it is convenient to express this system using matrix notation
[@Sanchez2013]. For this we define $\mathbf{P}(m, p) = (P_A(m, p), P_I(m, p),
P_R(m, p))^T$. Then the system of ODEs can be expressed as
\begin{equation}
\begin{split}
    \frac{d \mathbf{P}(m, p)}{dt} &= \mathbf{K} \mathbf{P}(m, p)
    - \mathbf{R}_m \mathbf{P}(m, p) + \mathbf{R}_m \mathbf{P}(m-1, p)
    - m \mathbf{\Gamma}_m \mathbf{P}(m, p) + (m + 1) \mathbf{\Gamma}_m \mathbf{P}(m + 1, p)\\
    &- m \mathbf{R}_p \mathbf{P}(m, p) + m \mathbf{R}_p \mathbf{P}(m, p - 1)
    - p \mathbf{\Gamma}_p \mathbf{P}(m, p) + (p + 1) \mathbf{\Gamma}_p \mathbf{P}(m, p + 1),
\end{split}
\end{equation}
where we defined matrices representing the promoter state transition
$\mathbf{K}$, 
\begin{equation}
 \mathbf{K} \equiv
 \begin{bmatrix}
    -k^{(p)}_{\text{off}}   & k^{(p)}_{\text{on}}         & 0\\
    k^{(p)}_{\text{off}}    & -k^{(p)}_{\text{on}} -k^{(r)}_{\text{on}}  & k^{(r)}_{\text{off}}\\
    0         & k^{(r)}_{\text{on}}         & -k^{(r)}_{\text{off}}
 \end{bmatrix},
\end{equation}
mRNA production, $\mathbf{R}_m$, and degradation, $\mathbf{\Gamma}_m$, as
\begin{equation}
\mathbf{R}_m \equiv
\begin{bmatrix}
    r_m   & 0 & 0\\
    0     & 0 & 0\\
    0     & 0 & 0\\
\end{bmatrix},
\end{equation}
and 
\begin{equation}
\mathbf{\Gamma}_m \equiv
  \begin{bmatrix}
    \gamma _m   & 0   & 0\\
    0     & \gamma _m & 0\\
    0     & 0   & \gamma _m\\
  \end{bmatrix}.
\end{equation}
For the protein we also define production $\mathbf{R}_p$ and degradation
$\mathbf{\Gamma}_p$ matrices as 
\begin{equation}
\mathbf{R}_p \equiv
\begin{bmatrix}
    r_p   & 0   & 0\\
    0     & r_p & 0\\
    0     & 0   & r_p\\
\end{bmatrix}
\end{equation}
and
\begin{equation}
\mathbf{\Gamma}_p \equiv
  \begin{bmatrix}
    \gamma _p   & 0   & 0\\
    0     & \gamma _p & 0\\
    0     & 0   & \gamma _p\\
\end{bmatrix}.
\end{equation}

The corresponding equation for the unregulated two-state promoter takes the
exact same form with the definition of the matrices following the same scheme
without including the third row and third column, and setting
$k^{(r)}_{\text{on}}$ and $k^{(r)}_{\text{off}}$ to zero.

A closed-form solution for this master equation might not even exist. The
approximate solution of chemical master equations of this kind is an active area
of research. As we will see in the two-state promoter master equation has been
analytically solved for the mRNA [@Peccoud1995] and protein distributions
[@Shahrezaei2008]. For our purposes, in we will detail how to use the Maximum
Entropy principle to approximate the full distribution for the two- and
three-state promoter.

## Parameter inference 

(Note: The Python code used for the calculations presented in this section can
be found in the [following
link](https://www.rpgroup.caltech.edu//chann_cap/software/chemical_master_mRNA_FISH_mcmc.html)
as an annotated Jupyter notebook)

With the objective of generating falsifiable predictions with meaningful
parameters, we infer the kinetic rates for this three-state promoter model using
different data sets generated in our lab over the last decade concerning
different aspects of the regulation of the simple repression motif. For example,
for the unregulated promoter transition rates $k^{(p)}_{\text{on}}$ and
$k^{(p)}_{\text{off}}$ and the mRNA production rate $r_m$, we use
single-molecule mRNA FISH counts from an unregulated promoter [@Jones2014a].
Once these parameters are fixed, we use the values to constrain the repressor
rates $k^{(r)}_{\text{on}}$ and $k^{(r)}_{\text{off}}$. These repressor rates
are obtained using information from mean gene expression measurements from bulk
LacZ colorimetric assays [@Garcia2011c]. We also expand our model to include the
allosteric nature of the repressor protein, taking advantage of video microscopy
measurements done in the context of multiple promoter copies [@Brewster2014] and
flow-cytometry measurements of the mean response of the system to different
levels of induction [@Razo-Mejia2018]. In what follows of this section we detail
the steps taken to infer the parameter values. At each step the values of the
parameters inferred in previous steps constrain the values of the parameters
that are not yet determined, building in this way a self-consistent model
informed by work that spans several experimental techniques.

### Unregulated promoter rates

We begin our parameter inference problem with the promoter on and off rates
$k^{(p)}_{\text{on}}$ and $k^{(p)}_{\text{off}}$, as well as the mRNA production
rate $r_m$. In this case there are only two states available to the promoter --
the inactive state $I$ and the transcriptionally active state $A$. That means
that the third ODE for $P_R(m, p)$ is removed from the system. The mRNA steady
state distribution for this particular two-state promoter model was solved
analytically by Peccoud and Ycart [@Peccoud1995]. This distribution $P(m) \equiv
P_I(m) + P_A(m)$ is of the form 
\begin{equation}
\small
  P(m \mid k^{(p)}_{\text{on}}, k^{(p)}_{\text{off}}, r_m, \gamma _m) =
  {\Gamma \left( \frac{k^{(p)}_{\text{on}}}{\gamma _m} + m \right) \over
  \Gamma(m + 1) \Gamma\left( \frac{k^{(p)}_{\text{off}}+k^{(p)}_{\text{on}}}{\gamma _m} + m \right)}
  {\Gamma\left( \frac{k^{(p)}_{\text{off}}+k^{(p)}_{\text{on}}}{\gamma _m} \right) \over
  \Gamma\left( \frac{k^{(p)}_{\text{on}}}{\gamma _m} \right)}
  \left( \frac{r_m}{\gamma _m} \right)^m
  F_1^1 \left( \frac{k^{(p)}_{\text{on}}}{\gamma _m} + m,
  \frac{k^{(p)}_{\text{off}} + k^{(p)}_{\text{on}}}{\gamma _m} + m,
  -\frac{r_m}{\gamma _m} \right),
\end{equation}
where $\Gamma(\cdot)$ is the gamma function, and $F_1^1$ is the confluent
hypergeometric function of the first kind. This rather complicated expression
will aid us to find parameter values for the rates. The inferred rates
$k^{(p)}_{\text{on}}$, $k^{(p)}_{\text{off}}$ and $r_m$ will be expressed in
units of the mRNA degradation rate $\gamma _m$. This is because the model in is
homogeneous in time, meaning that if we divide all rates by a constant it would
be equivalent to multiplying the characteristic time scale by the same constant.
As we will discuss in the next section, has degeneracy in the parameter values.
What this means is that a change in one of the parameters, specifically $r_m$,
can be compensated by a change in another parameter, specifically
$k^{(p)}_{\text{off}}$, to obtain the exact same distribution. To work around
this intrinsic limitation of the model we will include in our inference prior
information from what we know from equilibrium-based models.

### Bayesian parameter inference of RNAP rates

In order to make progress at inferring the unregulated promoter state transition
rates, we make use of the single-molecule mRNA FISH data from Jones et al.
[@Jones2014a]. shows the distribution of mRNA per cell for the *lacUV5* promoter
used for our inference. This promoter, being very strong, has a mean copy number
of $\langle m \rangle \approx 18$ mRNA/cell.

![***lacUV5* mRNA per cell distribution.** Data from [@Jones2014a] of the
unregulated *lacUV5* promoter as inferred from single molecule mRNA
FISH.](ch5_fig01){#fig:ch5_fig01 short-caption="lacUV5* mRNA per cell
distribution"}

Having this data in hand we now turn to Bayesian parameter inference. Writing
Bayes theorem we have
\begin{equation}
P(k^{(p)}_{\text{on}}, k^{(p)}_{\text{off}}, r_m \mid D) = 
\frac{P(D \mid k^{(p)}_{\text{on}}, k^{(p)}_{\text{off}}, r_m)
P(k^{(p)}_{\text{on}}, k^{(p)}_{\text{off}}, r_m)}{P(D)},
\end{equation}
where $D$ represents the data. For this case the data consists of single-cell
mRNA counts $D = \{ m_1, m_2, \ldots, m_N \}$, where $N$ is the number of cells.
We assume that each cell's measurement is independent of the others such that we
can rewrite as 
\begin{equation}
P(k^{(p)}_{\text{on}}, k^{(p)}_{\text{off}}, r_m \mid \{m_i\}) \propto
  \left[\prod_{i=1}^N P(m_i \mid k^{(p)}_{\text{on}}, k^{(p)}_{\text{off}}, r_m) \right]
  P(k^{(p)}_{\text{on}}, k^{(p)}_{\text{off}}, r_m),
\end{equation}
where we ignore the normalization constant $P(D)$. The likelihood term $P(m_i
\mid k^{(p)}_{\text{on}}, k^{(p)}_{\text{off}}, r_m)$ is exactly given by with
$\gamma _m = 1$. Given that we have this functional form for the distribution,
we can use Markov Chain Monte Carlo (MCMC) sampling to explore the 3D parameter
space in order to fit to the mRNA-FISH data.

#### Constraining the rates given prior thermodynamic knowledge. 

One of the strengths of the Bayesian approach is that we can include all the
prior knowledge on the parameters when performing an inference [@MacKay2003].
Basic features such as the fact that the rates have to be strictly positive
constrain the values that these parameters can take. For the specific rates
analyzed in this section we know more than the simple constraint of non-negative
values. The expression of an unregulated promoter has been studied from a
thermodynamic perspective [@Brewster2012]. Given the underlying assumptions of
these equilibrium models, in which the probability of finding the RNAP bound to
the promoter is proportional to the transcription rate [@Bintu2005a], they can
only make statements about the mean expression level. Nevertheless if both the
thermodynamic and the kinetic model describe the same process, the predictions
for the mean gene expression level must agree. That means that we can use what
we know about the mean gene expression, and how this is related to parameters
such as molecule copy numbers and binding affinities, to constrain the values
that the rates in question can take.

In the case of this two-state promoter it can be shown that the mean number of
mRNA is given by [@Sanchez2013] (See XXX for moment computation)
\begin{equation}
\langle m \rangle = \frac{r_m}{\gamma _m}
\frac{k^{(p)}_{\text{on}}}{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}.
\end{equation}
Another way of expressing this is as $\frac{r_m}{\gamma _m} \times
p_{\text{active}}^{(p)}$, where $p_{\text{active}}^{(p)}$ is the probability of
the promoter being in the transcriptionally active state. The thermodynamic
picture has an equivalent result where the mean number of mRNA is given by
[@Brewster2012; @Bintu2005a]
\begin{equation}
\left\langle m \right\rangle = \frac{r_m}{\gamma _m}
\frac{\frac{P}{N_{NS}} e^{-\beta\Delta\varepsilon_p}
}{1 + \frac{P}{N_{NS}} e^{-\beta\Delta\varepsilon_p}},
\end{equation}
where $P$ is the number of RNAP per cell, $N_{NS}$ is the number of non-specific
binding sites, $\Delta\varepsilon_p$ is the RNAP binding energy in $k_BT$ units
and $\beta\equiv {(k_BT)}^{-1}$. Using and we can easily see that if these
frameworks are to be equivalent, then it must be true that
\begin{equation}
\frac{k^{(p)}_{\text{on}}}{ k^{(p)}_{\text{off}}} =
\frac{P}{N_{NS}} e^{-\beta\Delta\varepsilon_p},
\end{equation}
or equivalently
\begin{equation}
\ln \left(\frac{k^{(p)}_{\text{on}}}{k^{(p)}_{\text{off}}}\right) =
  -\beta\Delta\varepsilon_p + \ln P - \ln N_{NS}.
\end{equation}
To put numerical values into these variables we can use information from the
literature. The RNAP copy number is order $P \approx 1000-3000$ RNAP/cell for a
1 hour doubling time [@Klumpp2008]. As for the number of non-specific binding
sites and the binding energy, we have that $N_{NS} = 4.6\times 10^6$
[@Bintu2005a] and $-\beta\Delta\varepsilon_p \approx 5 - 7$ [@Brewster2012].
Given these values we define a Gaussian prior for the log ratio of these two
quantities of the form
\begin{equation}
P\left(\ln \left(\frac{k^{(p)}_{\text{on}}}{k^{(p)}_{\text{off}}}\right) \right) \propto
\exp \left\{ - \frac{\left(\ln \left(\frac{k^{(p)}_{\text{on}}}{k^{(p)}_{\text{off}}}\right) -
\left(-\beta\Delta\varepsilon_p + \ln P - \ln N_{NS} \right) \right)^2
}{2 \sigma^2} \right\},
\end{equation}
where $\sigma$ is the variance that accounts for the uncertainty in these
parameters. We include this prior as part of the prior term
$P(k^{(p)}_{\text{on}}, k^{(p)}_{\text{off}}, r_m)$ of . We then use MCMC to
sample out of the posterior distribution given by . shows the MCMC samples of
the posterior distribution. For the case of the $k^{(p)}_{\text{on}}$ parameter
there is a single symmetric peak. $k^{(p)}_{\text{off}}$ and $r_m$ have a rather
long tail towards large values. In fact, the 2D projection of
$k^{(p)}_{\text{off}}$ vs $r_m$ shows that the model is sloppy, meaning that the
two parameters are highly correlated. This feature is a common problem for many
non-linear systems used in biophysics and systems biology [@Transtrum2015]. What
this implies is that we can change the value of $k^{(p)}_{\text{off}}$, and then
compensate by a change in $r_m$ in order to maintain the shape of the mRNA
distribution. Therefore it is impossible from the data and the model themselves
to narrow down a single value for the parameters. Nevertheless since we included
the prior information on the rates as given by the analogous form between the
equilibrium and non-equilibrium expressions for the mean mRNA level, we obtained
a more constrained parameter value for the RNAP rates and the transcription rate
that we will take as the peak of this long-tailed distribution.

![**MCMC posterior distribution.** Sampling out of the plot shows 2D and 1D
projections of the 3D parameter space. The parameter values are (in units of the
mRNA degradation rate $\gamma _m$) $k^{(p)}_{\text{on}} = 4.3^{+1}_{-0.3}$,
$k^{(p)}_{\text{off}} = 18.8^{+120}_{-10}$ and $r_m = 103.8^{+423}_{-37}$ which
are the modes of their respective distributions, where the superscripts and
subscripts represent the upper and lower bounds of the 95$^\text{th}$ percentile
of the parameter value distributions](ch5_fig02){#fig:ch5_fig02
short-caption="MCMC posterior distribution."}

The inferred values $k^{(p)}_{\text{on}} = 4.3^{+1}_{-0.3}$,
$k^{(p)}_{\text{off}} = 18.8^{+120}_{-10}$ and $r_m = 103.8^{+423}_{-37}$ are
given in units of the mRNA degradation rate $\gamma _m$. Given the asymmetry of
the parameter distributions we report the upper and lower bound of the
95$^\text{th}$ percentile of the posterior distributions. Assuming a mean
life-time for mRNA of $\approx$ 3 min (from this
[link](http://bionumbers.hms.harvard.edu/bionumber.aspx?&id=107514&ver=1&trm=mRNA%20mean%20lifetime))
we have an mRNA degradation rate of $\gamma _m \approx 2.84 \times 10^{-3}
s^{-1}$. Using this value gives the following values for the inferred rates:
$k^{(p)}_{\text{on}} = 0.024_{-0.002}^{+0.005} s^{-1}$, $k^{(p)}_{\text{off}} =
{0.11}_ {-0.05}^{+0.66} s^{-1}$, and $r_m = 0.3_{-0.2}^{+2.3} s^{-1}$.

[@Fig:ch5_fig03] compares the experimental data from with the resulting
distribution obtained by substituting the most likely parameter values into . As
we can see this two-state model fits the data adequately.

![**Experimental vs. theoretical distribution of mRNA per cell using parameters
from Bayesian inference.** Dotted line shows the result of using along with the
parameters inferred for the rates. Blue bars are the same data as obtained from
[@Jones2014a].](ch5_fig03){#fig:ch5_fig03 short-caption="Experimental vs.
theoretical distribution of mRNA per cell using parameters from Bayesian
inference"}

### Accounting for variability in the number of promoters

As discussed in ref. [@Jones2014a] and further expanded in [@Peterson2015] an
important source of cell-to-cell variability in gene expression in bacteria is
the fact that, depending on the growth rate and the position relative to the
chromosome replication origin, cells can have multiple copies of any given gene.
Genes closer to the replication origin have on average higher gene copy number
compared to genes at the opposite end. For the locus in which our reporter
construct is located (*galK*) and the doubling time of the mRNA FISH experiments
we expect to have $\approx$ 1.66 copies of the gene [@Jones2014a; @Bremer1996].
This implies that the cells spend 2/3 of the cell cycle with two copies of the
promoter and the rest with a single copy.

To account for this variability in gene copy we extend the model assuming that
when cells have two copies of the promoter the mRNA production rate is $2 r_m$
compared to the rate $r_m$ for a single promoter copy. The probability of
observing a certain mRNA copy $m$ is therefore given by
\begin{equation}
P(m) = P(m \mid \text{one promoter}) \cdot P(\text{one promoter}) +
P(m \mid \text{two promoters}) \cdot P(\text{two promoters}).
\end{equation}
Both terms $P(m \mid \text{promoter copy})$ are given by with the only
difference being the rate $r_m$. It is important to acknowledge that assumes
that once the gene is replicated the time scale in which the mRNA count relaxes
to the new steady state is much shorter than the time that the cells spend in
this two promoter copies state. This approximation should be valid for a short
lived mRNA molecule, but the assumption is not applicable for proteins whose
degradation rate is comparable to the cell cycle length as explored in XXX.

In order to repeat the Bayesian inference including this variability in gene
copy number we must split the mRNA count data into two sets -- cells with a
single copy of the promoter and cells with two copies of the promoter. For the
single molecule mRNA FISH data there is no labeling of the locus, making it
impossible to determine the number of copies of the promoter for any given cell.
We therefore follow Jones et al. [@Jones2014a] in using the cell area as a proxy
for stage in the cell cycle. In their approach they sorted cells by area,
considering cells below the 33$^\text{th}$ percentile having a single promoter
copy and the rest as having two copies. This approach ignores that cells are not
uniformly distributed along the cell cycle. As first derived in [@Powell1956]
populations of cells in a log-phase are exponentially
distributed along the cell cycle. This distribution is of the form
\begin{equation}
P(a) = (\ln 2) \cdot 2^{1 - a},
\end{equation}
where $a \in [0, 1]$ is the stage of the cell cycle, with $a = 0$ being the
start of the cycle and $a = 1$ being the cell division (See for a derivation of
). shows the separation of the two groups based on area where was used to weight
the distribution along the cell cycle.

![**Separation of cells based on cell size.** Using the area as a proxy for
position in the cell cycle, cells can be sorted into two groups -- small cells
(with one promoter copy) and large cells (with two promoter copies). The
vertical black line delimits the threshold that divides both groups as weighted
by .](ch5_fig04){#fig:ch5_fig04 short-caption="Separation of cells based on cell
size"}

A subtle, but important consequence of is that computing any quantity for a
single cell is not equivalent to computing the same quantity for a population of
cells. For example, let us assume that we want to compute the mean mRNA copy
number $\langle m \rangle$. For a single cell this would be of the form
\begin{equation}
\langle m \rangle_{\text{cell}} =
\langle m \rangle_1 \cdot f + \langle m \rangle_2 \cdot (1 - f),
\end{equation}
where $\langle m \rangle_i$ is the mean mRNA copy number with $i$ promoter
copies in the cell, and $f$ is the fraction of the cell cycle that cells spend
with a single copy of the promoter. For a single cell the probability of having
a single promoter copy is equivalent to this fraction $f$. But tells us that if
we sample unsynchronized cells we are not sampling uniformly across the cell
cycle. Therefore for a population of cells the mean mRNA is given by
\begin{equation}
\langle m \rangle_{\text{population}} =
\langle m \rangle_1 \cdot \phi + \langle m \rangle_2 \cdot (1 - \phi)
\end{equation}
where the probability of sampling a cell with one promoter $\phi$ is given by
\begin{equation}
\phi = \int_0^f P(a) da,
\end{equation}
where $P(a)$ is given by XXX. What this equation computes is the probability of
sampling a cell during a stage of the cell cycle $< f$ where the reporter gene
hasn't been replicated yet. shows the distribution of both groups. As expected
larger cells have a higher mRNA copy number on average.

![**mRNA distribution for small and large cells.** (A) histogram and (B)
cumulative distribution function of the small and large cells as determined in .
The triangles above histograms in (A) indicate the mean mRNA copy number for
each group.](ch5_fig05){#fig:ch5_fig05 short-caption="mRNA distribution for
small and large cells."}

We modify to account for the two separate groups of cells. Let $N_s$ be the
number of cells in the small size group and $N_l$ the number of cells in the
large size group. Then the posterior distribution for the parameters is of the
form 
\begin{equation}
\small
P(k^{(p)}_{\text{on}}, k^{(p)}_{\text{off}}, r_m \mid \{m_i\}) \propto
\left[\prod_{i=1}^{N_s} P(m_i \mid k^{(p)}_{\text{on}}, k^{(p)}_{\text{off}}, r_m)\right]
\left[\prod_{j=1}^{N_l} P(m_j \mid k^{(p)}_{\text{on}}, k^{(p)}_{\text{off}}, 2 r_m)\right]
P(k^{(p)}_{\text{on}}, k^{(p)}_{\text{off}}, r_m),
\end{equation}
where we split the product of small and large cells.

For the two-promoter model the prior shown in requires a small modification.
gives the mean mRNA copy number of a population of asynchronous cells growing at
steady state. Given that we assume that the only difference between having one
vs. two promoter copies is the change in transcription rate from $r_m$ in the
single promoter case to $2 r_m$ in the two-promoter case we can write as
\begin{equation}
\langle m \rangle = \phi \cdot \frac{r_m}{\gamma _m}
\frac{k^{(p)}_{\text{on}}}{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}} +
(1 -\phi) \cdot \frac{2 r_m}{\gamma _m}
\frac{k^{(p)}_{\text{on}}}{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}.
\end{equation}
This can be simplified to
\begin{equation}
\langle m \rangle = (2 - \phi) \frac{r_m}{\gamma _m} 
\frac{k^{(p)}_{\text{on}}}{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}.
\end{equation}

Equating and to again require self-consistent predictions of the mean mRNA from
the equilibrium and kinetic models gives
\begin{equation}
(2 - \phi) \frac{k^{(p)}_{\text{on}}}{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}} =
\frac{\frac{P}{N_{NS}} e^{-\beta\Delta\varepsilon_p}
}{1 + \frac{P}{N_{NS}} e^{-\beta\Delta\varepsilon_p}}.
\end{equation}
Solving for
$\frac{k^{(p)}_{\text{on}}}{k^{(p)}_{\text{off}}}$ results in 
\begin{equation}
\left(\frac{k^{(p)}_{\text{on}}}{k^{(p)}_{\text{off}}}\right) =
\frac{\rho}{\left[ (1 + \rho)(2 - \phi) - \rho \right]},
\end{equation}
where we define $\rho \equiv \frac{P}{N_{NS}} e^{-\beta\Delta\varepsilon_p}$. To
simplify things further we notice that for the specified values of $P = 1000 -
3000$ per cell, $N_{NS} = 4.6 \times 10^6$ bp, and $-\beta\Delta\varepsilon_p =
5 - 7$, we can safely assume that $\rho \ll 1$. This simplifying assumption has
been previously called the weak promoter approximation [@Garcia2011c]. Given
this we can simplify as
\begin{equation}
\frac{k^{(p)}_{\text{on}}}{k^{(p)}_{\text{off}}} =
\frac{1}{2 - \phi} \frac{P}{N_{NS}} e^{-\beta\Delta\varepsilon_p}.
\end{equation}
Taking the log of both sides gives
\begin{equation}
\ln\left(\frac{k^{(p)}_{\text{on}}}{k^{(p)}_{\text{off}}}\right) =
-\ln (2 - \phi) + \ln P - \ln N_{NS}
  - \beta\Delta\varepsilon_p.
\end{equation}
With this we can set as before a Gaussian prior to
constrain the ratio of the RNAP rates as
\begin{equation}
P\left(\ln \left(\frac{k^{(p)}_{\text{on}}}{k^{(p)}_{\text{off}}}\right) \right)  \propto
\exp \left\{ - \frac{\left(\ln \left(\frac{k^{(p)}_{\text{on}}}{k^{(p)}_{\text{off}}}\right) -
\left[-\ln(2 - \phi) -\beta\Delta\varepsilon_p + \ln P - \ln N_{NS} \right) \right]^2
}{2 \sigma^2} \right\}.
\end{equation}

[@Fig:ch5_fig06] shows the result of sampling out of XXX. Again we see that the
model is highly sloppy with large credible regions obtained for
$k^{(p)}_{\text{off}}$ and $r_m$. Nevertheless, again the use of the prior
information allows us to get a parameter values consistent with the equilibrium
picture.

![**MCMC posterior distribution for a multi-promoter model.** Sampling out of
the plot shows 2D and 1D projections of the 3D parameter space. The parameter
values are (in units of the mRNA degradation rate $\gamma _m$)
$k^{(p)}_{\text{on}} = 6.4^{+0.8}_{-0.4}$, $k^{(p)}_{\text{off}} =
132^{+737}_{-75}$ and $r_m = 257^{+1307}_{-132}$ which are the modes of their
respective distributions, where the superscripts and subscripts represent the
upper and lower bounds of the 95$^\text{th}$ percentile of the parameter value
distributions. The sampling was bounded to values $<$ 1000 for numerical
stability when computing the confluent hypergeometric
function.](ch5_fig06){#fig:ch5_fig06 short-caption="MCMC posterior distribution
for a multi-promoter model"}

Using again a mRNA mean lifetime of $\approx 3$ min gives the following values
for the parameters: $k^{(p)}_{\text{on}} = {0.03}_{-0.002}^{+0.004} s^{-1}$,
$k^{(p)}_{\text{off}} = {0.7} _{-0.4}^{+4.1} s^{-1}$, and $r_m =
{1.4}_{-0.7}^{+7.3} s^{-1}$. shows the result of applying using these parameter
values. Specifically (A) shows the global distribution including cells with one
and two promoters and (B) splits the distributions within the two populations.
Given that the model adequately describes both populations independently and
pooled together we confirm that using the cell area as a proxy for stage in the
cell cycle and the doubling of the transcription rate once cells have two
promoters are reasonable approximations.

![**Experimental vs. theoretical distribution of mRNA per cell using parameters
for multi-promoter model.** (A) Solid line shows the result of using with the
parameters inferred by sampling . Blue bars are the same data as from
[@Jones2014a]. (B) Split distributions of small cells (light blue bars) and
large cells (dark blue) with the corresponding theoretical predictions with
transcription rate $r_m$ (light blue line) and transcription rate $2 r_m$ (dark
blue line)](ch5_fig07){#fig:ch5_fig07 short-caption="Experimental vs.
theoretical distribution of mRNA per cell using parameters for multi-promoter
model"}

It is hard to make comparisons with literature reported values because these
kinetic rates are effective parameters hiding a lot of the complexity of
transcription initiation [@Browning2004]. Also the non-identifiability of the
parameters restricts our explicit comparison of the actual numerical values of
the inferred rates. Nevertheless from the model we can see that the mean burst
size for each transcription event is given by $r_m / k^{(p)}_{\text{off}}$. From
our inferred values we obtain then a mean burst size of $\approx 1.9$
transcripts per cell. This is similar to the reported burst size of 1.15 on a
similar system on *E. coli* [@Yu2006].

### Repressor rates from three-state regulated promoter.

Having determined the unregulated promoter transition rates we now proceed to
determine the repressor rates $k^{(r)}_{\text{on}}$ and $k^{(r)}_{\text{off}}$.
The values of these rates are constrained again by the correspondence between
our kinetic picture and what we know from equilibrium models [@Phillips2015a].
For this analysis we again exploit the feature that, at the mean, both the
kinetic language and the thermodynamic language should have equivalent
predictions. Over the last decade there has been great effort in developing
equilibrium models for gene expression regulation [@Buchler2003; @Vilar2011;
@Bintu2005a]. In particular our group has extensively characterized the simple
repression motif using this formalism [@Garcia2011c; @Brewster2014;
@Razo-Mejia2018].

The dialogue between theory and experiments has led to simplified expressions
that capture the phenomenology of the gene expression response as a function of
natural variables such as molecule count and affinities between molecular
players. A particularly interesting quantity for the simple repression motif
used by Garcia & Phillips [@Garcia2011c] is the fold-change in gene expression,
defined as
\begin{equation}
\text{fold-change} =
\frac{\left\langle{\text{gene expression}(R \neq 0)}\right\rangle}{
\left\langle{\text{gene expression}(R = 0)}\right\rangle},
\end{equation}
where $R$ is the number of repressors per cell and
$\left\langle{\cdot}\right\rangle$ is the population average. The fold-change is
simply the mean expression level in the presence of the repressor relative to
the mean expression level in the absence of regulation. In the language of
statistical mechanics this quantity takes the form
\begin{equation}
\text{fold-change} = \left( 1 + \frac{R}{N_{NS}}
e^{-\beta\Delta\varepsilon_r} \right)^{-1},
\end{equation}
where $\Delta\varepsilon_r$ is the repressor-DNA binding energy, and as before
$N_{NS}$ is the number of non-specific binding sites where the repressor can
bind [@Garcia2011c].

To compute the fold-change in the chemical master equation language we compute
the first moment of the steady sate mRNA distribution $\langle m \rangle$ for
both the three-state promoter ($R \neq 0$) and the two-state promoter case
($R=0$) (See for moment derivation). The unregulated (two-state) promoter mean
mRNA copy number is given by . For the regulated (three-state) promoter we have
an equivalent expression of the form 
\begin{equation}
\left\langle{m (R \neq 0)}\right\rangle = 
(2 - \phi)\frac{r_m}{\gamma _m}
\frac{k^{(r)}_{\text{off}}k^{(p)}_{\text{on}}
}{k^{(p)}_{\text{off}}k^{(r)}_{\text{off}} +
k^{(p)}_{\text{off}}k^{(r)}_{\text{on}} + 
k^{(r)}_{\text{off}}k^{(p)}_{\text{on}}}.
\end{equation}
Computing the fold-change then gives
\begin{equation}
\text{fold-change} =
\frac{\left\langle{m (R \neq 0)}\right\rangle}
{\left\langle{m (R = 0)}\right\rangle} =
\frac{k^{(r)}_{\text{off}} \left( k^{(p)}_{\text{off}} + 
k^{(p)}_{\text{on}} \right)}{
k^{(p)}_{\text{off}}k^{(r)}_{\text{on}} +
k^{(r)}_{\text{off}} \left( k^{(p)}_{\text{off}} + k^{(p)}_{\text{on}} \right)},
\end{equation}
where the factor $(2 - \phi)$ due to the multiple promoter copies, the
transcription rate $r_m$ and the mRNA degradation rate $\gamma _m$ cancel out.

Given that the number of repressors per cell $R$ is an experimental variable
that we can control, we assume that the rate at which the promoter transitions
from the transcriptionally inactive state to the repressor bound state,
$k^{(r)}_{\text{on}}$, is given by the concentration of repressors $[R]$ times a
diffusion limited on rate $k_o$. For the diffusion limited constant $k_o$ we use
the value used by Jones et al. [@Jones2014a]. With this in hand we can rewrite
as
\begin{equation}
\text{fold-change} = \left( 1 + \frac{k_0 [R]}{k^{(r)}_{\text{off}}}
\frac{k^{(p)}_{\text{off}}}{k^{(p)}_{\text{on}} +
k^{(p)}_{\text{off}}} \right)^{-1}.
\end{equation}

We note that both and have the same functional form. Therefore if both languages
predict the same output for the mean gene expression level, it must be true that
\begin{equation}
\frac{k_o [R]}{k^{(r)}_{\text{off}}}
\frac{k^{(p)}_{\text{off}}}{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}} =
\frac{R}{N_{NS}} e^{-\beta\Delta\varepsilon_r}.
\end{equation}
Solving for $k^{(r)}_{\text{off}}$ gives
\begin{equation}
k^{(r)}_{\text{off}} = \frac{k_o [R] N_{NS}
e^{\beta\Delta\varepsilon_r}}{R}
\frac{k^{(p)}_{\text{off}}}{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}.
\end{equation}

Since the reported value of $k_o$ is given in units of nM$^{-1}$s$^{-1}$
in order for the units to cancel properly the repressor concentration
has to be given in nM rather than absolute count. If we consider that
the repressor concentration is equal to
\begin{equation}
[R] = \frac{R}{V_{cell}}\cdot \frac{1}{N_A},
\end{equation}
where $R$ is the absolute repressor copy number per cell, $V_{cell}$ is the cell
volume and $N_A$ is Avogadro's number. The *E. coli* cell volume is 2.1 fL
[@Radzikowski2016], and Avogadro's number is $6.022 \times 10^{23}$. If we
further include the conversion factor to turn M into nM we find that
\begin{equation}
[R] = \frac{R}{2.1 \times 10^{-15} L} \cdot \frac{1}{6.022 \times 10^{23}}
\cdot \frac{10^9 \text{ nmol}}{1 \text{ mol}} \approx 0.8 \times R.
\end{equation}
Using this we simplify as
\begin{equation}
k^{(r)}_{\text{off}} \approx 0.8 \cdot k_o \cdot
N_{NS} e^{\beta\Delta\varepsilon_r}
\cdot \frac{k^{(p)}_{\text{off}}}{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}.
\end{equation}
What shows is the direct relationship that must be satisfied if the equilibrium
model is set to be consistent with the non-equilibrium kinetic picture.
summarizes the values obtained for the three operator sequences used throughout
this work. To compute these numbers the number of non-specific binding sites
$N_{NS}$ was taken to be $4.6 \times 10^6$ bp, i.e. the size of the *E. coli*
K12 genome.

| Operator | $\Delta\varepsilon_r\; (k_BT)$ | $k^{(r)}_{\text{off}}\; (s^{-1})$ |
| -------- | ------------------------------ | --------------------------------- |
| O1       | -15.3                          | 0.002                             |
| O2       | -13.9                          | 0.008                             |
| O3       | -9.7                           | 0.55                              |
Table: **Binding sites and corresponding parameters.** {#tbl:ch5_tbl01}


*In-vivo* measurements of the Lac repressor off rate have been done with
single-molecule resolution [@Hammar2014]. The authors report a mean residence
time of $5.3 \pm 0.2$ minutes for the repressor on an O1 operator. The
corresponding rate is $k^{(r)}_{\text{off}} \approx 0.003$ $(s^{-1})$, very
similar value to what we inferred from our model. In this same reference the
authors determined that on average the repressor takes $30.9 \pm 0.5$ seconds to
bind to the operator [@Hammar2014]. Given the kinetic model presented in (A)
this time can be converted to the probability of not being on the repressor
bound state $P_{\text{not }R}$. This is computed as
\begin{equation}
P_{\text{not }R} = {\tau_{\text{not }R} \over
\tau_{\text{not }R} + \tau_{R}},
\end{equation}
where $\tau_{\text{not }R}$ is the average time that the operator is not
occupied by the repressor and $\tau_{R}$ is the average time that the repressor
spends bound to the operator. Substituting the numbers from [@Hammar2014] gives
$P_{\text{not }R} \approx 0.088$. From our model we can compute the zeroth
moment $\left\langle{m^0 p^0}\right\rangle$ for each of the three promoter
states. This moment is equivalent to the probability of being on each of the
promoter states. Upon substitution of our inferred rate parameters we can
compute $P_{\text{not }R}$ as 
\begin{equation}
P_{\text{not }R} = 1 - P_R \approx 0.046,
\end{equation}
where $P_R$ is the probability of the promoter being bound by the repressor. The
value we obtained is within a factor of two from the one reported in
[@Hammar2014].
