## Results

### Minimal model of transcriptional regulation {#sec_model}

As a tractable circuit for which we have control over the parameters both
theoretically and experimentally, we chose the so-called simple repression
motif, a common regulatory scheme among prokaryotes [@Rydenfelt2014]. This
circuit consists of a single promoter with an RNA-polymerase (RNAP) binding site
and a single binding site for a transcriptional repressor [@Garcia2011c]. The
regulation due to the repressor occurs via exclusion of the RNAP from its
binding site when the repressor is bound, decreasing the likelihood of having a
transcription event. As with many important macromolecules, we consider the
repressor to be allosteric, meaning that it can exist in two conformations, one
in which the repressor is able to bind to the specific binding site (active
state) and one in which it cannot bind the specific binding site (inactive
state). The environmental signaling occurs via passive import of an
extracellular inducer that binds the repressor, shifting the equilibrium between
the two conformations of the repressor [@Razo-Mejia2018]. In previous work we
have extensively characterized the mean response of this circuit under different
conditions using equilibrium based models [@Phillips2019]. Here we build upon
these models to characterize the full distribution of gene expression with
parameters such as repressor copy number and its affinity for the DNA being
systematically varied.

As the copy number of molecular species is a discrete quantity, chemical master
equations have emerged as a useful tool to model their inherent probability
distribution [@Sanchez2013]. In [@Fig:ch3_fig02](A) we show the minimal model
and the necessary set of parameters needed to compute the full distribution of
mRNA and its protein gene product. Specifically, we assume a three-state model
where the promoter can be found in a 1) transcriptionally active state ($A$
state), 2) a transcriptionally inactive state without the repressor bound ($I$
state) and 3) a transcriptionally inactive state with the repressor bound ($R$
state). We do not assume that the transition between the active state $A$ and
the inactive state $I$ occurs due to RNAP binding to the promoter as the
transcription initiation kinetics involve several more steps than simple binding
[@Browning2004]. We coarse-grain all these steps into effective "*on*" and
"*off*" states for the promoter, consistent with experiments demonstrating the
bursty nature of gene expression in *E. coli* [@Golding2005]. These three states
generate a system of coupled differential equations for each of the three state
distributions $P_A(m, p; t)$, $P_I(m, p; t)$ and $P_R(m, p; t)$, where $m$ and
$p$ are the mRNA and protein count per cell, respectively and $t$ is time. Given
the rates depicted in [@Fig:ch3_fig02](A) we define the system of ODEs for a
specific $m$ and $p$. For the transcriptionally active state, we have
$$
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
    - \overbrace{\gamma _p p P_A(m, p)}^{p \rightarrow p-1}, % p -> p-1
\end{split}
$${#eq:ch3_eq02}
where the state transitions for each term are labeled by overbraces. For the
transcriptionally inactive state $I$, we have
$$
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
$${#eq:ch3_eq03}
And finally, for the repressor bound state $R$,
$$
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
$${#eq:ch3_eq04}
As we will discuss later, the protein degradation term $\gamma _p$ is set to
zero since active protein degradation is slow compared to the cell cycle of
exponentially growing bacteria, but rather we explicitly implement binomial
partitioning of the proteins into daughter cells upon division [@Maurizi1992].

It is convenient to rewrite these equations in a compact matrix notation
[@Sanchez2013]. For this we define the vector $\mathbf{P}(m, p)$ as 
$$
\mathbf{P}(m, p) = (P_A(m, p), P_I(m, p), P_R(m, p))^T,
$${#eq:ch3_eq05}
where $^T$ is the transpose. By defining the matrices $\mathbf{K}$ to contain
the promoter state transitions, $\mathbf{R}_m$ and $\mathbf{\Gamma} _m$ to
contain the mRNA production and degradation terms, respectively, and
$\mathbf{R}_p$ and $\mathbf{\Gamma}_p$ to contain the protein production and
degradation terms, respectively, the system of ODEs can then be written as (See
for full definition of these matrices)
$$
\begin{split}
    \frac{\mathbf{P}(m, p)}{dt} &= 
    \left( \mathbf{K} - \mathbf{R}_m -m\mathbf{\Gamma} _m 
    -m \mathbf{R}_p -p\mathbf{\Gamma}_p \right) \mathbf{P}(m, p)\\
    &+ \mathbf{R}_m \mathbf{P}(m-1, p)
    + (m + 1) \mathbf{\Gamma} _m \mathbf{P}(m + 1, p)\\
    &+ m \mathbf{R}_p \mathbf{P}(m, p - 1)
    + (p + 1) \mathbf{\Gamma}_p \mathbf{P}(m, p + 1).
\end{split}
$${#eq:ch3_eq06}
Having defined the gene expression dynamics we now proceed to determine all rate
parameters in [@Eq:ch3_eq06].

### Inferring parameters from published data sets {#sec_param}

A decade of research in our group has characterized the simple repression motif
with an ever expanding array of predictions and corresponding experiments to
uncover the physics of this genetic circuit [@Phillips2019]. In doing so we have
come to understand the mean response of a single promoter in the presence of
varying levels of repressor copy numbers and repressor-DNA affinities
[@Garcia2011c], due to the effect that competing binding sites and multiple
promoter copies impose [@Brewster2014], and in recent work, assisted by the
Monod-Wyman-Changeux (MWC) model, we expanded the scope to the allosteric nature
of the repressor [@Razo-Mejia2018]. All of these studies have exploited the
simplicity and predictive power of equilibrium approximations to these
non-equilibrium systems [@Buchler2003]. We have also used a similar kinetic
model to that depicted in [@Fig:ch3_fig02](A) to study the noise in mRNA copy
number [@Jones2014a]. Although these studies focus on the same experimental
system described by different theoretical frameworks, in earlier work in our
laboratory an attempt to unite parametric knowledge across studies based on
equilibrium and non-equilibrium models has not been performed previously. As a
test case of the depth of our theoretical understanding of this simple
transcriptional regulation system we combine all of the studies mentioned above
to inform the parameter values of the model presented in [@Fig:ch3_fig02](A).
[@Fig:ch3_fig02](B) schematizes the data sets and experimental techniques used
to measure gene expression along with the parameters that can be inferred from
them.

Section XXX expands on the details of how the inference was performed for each
of the parameters. Briefly, the promoter activation and inactivation rates
$k^{(p)}_{\text{on}}$ and $k^{(p)}_{\text{off}}$, as well as the transcription
rate $r_m$ were obtained in units of the mRNA degradation rate $\gamma _m$ by
fitting a two-state promoter model (no state $R$ from [@Fig:ch3_fig02](A))
[@Peccoud1995] to mRNA FISH data of an unregulated promoter (no repressor
present in the cell) [@Jones2014a]. The repressor on rate is assumed to be of
the form $k^{(r)}_{\text{on}} = k_o [R]$ where $k_o$ is a diffusion-limited on
rate and $[R]$ is the concentration of active repressor in the cell
[@Jones2014a]. This concentration of active repressor is at the same time
determined by the repressor copy number in the cell, and the fraction of these
repressors that are in the active state, i.e. able to bind DNA. Existing
estimates of the transition rates between conformations of allosteric molecules
set them at the microsecond scale [@Cui2008]. By considering this to be
representative for our repressor of interest, the separation of time-scales
between the rapid conformational changes of the repressor and the slower
downstream processes such as the open-complex formation processes allow us to
model the probability of the repressor being in the active state as an
equilibrium MWC process. The parameters of the MWC model $K_A$, $K_I$ and
$\Delta\varepsilon_{AI}$ were previously characterized from video-microscopy and
flow-cytometry data [@Razo-Mejia2018]. For the repressor off rate,
$k^{(r)}_{\text{off}}$, we take advantage of the fact that the mean mRNA copy
number as derived from the model in [@Fig:ch3_fig02](A) cast in the language of
rates is of the same functional form as the equilibrium model cast in the
language of binding energies [@Phillips2015a]. Therefore the value of the
repressor-DNA binding energy $\Delta\varepsilon_r$ constrains the value of the
repressor off rate $k^{(r)}_{\text{off}}$. These constraints on the rates allow
us to make self-consistent predictions under both the equilibrium and the
kinetic framework. Having all parameters in hand, we can now proceed to solve
the gene expression dynamics.

![**Minimal kinetic model of transcriptional regulation for a simple repression
architecture.** (A) Three-state promoter stochastic model of transcriptional
regulation by a repressor. The regulation by the repressor occurs via exclusion
of the transcription initiation machinery, not allowing the promoter to
transition to the transcriptionally active state. All parameters highlighted
with colored boxes were determined from published datasets based on the same
genetic circuit. Parameters in dashed boxes were taken directly from values
reported in the literature or adjusted to satisfy known biological restrictions.
(B) Data sets used to infer the parameter values. From left to right Garcia &
Phillips [@Garcia2011c] is used to determine $k^{(r)}_{\text{off}}$ and
$k^{(r)}_{\text{on}}$, Brewster et al. [@Brewster2014] is used to determine
$\Delta\varepsilon_{AI}$ and $k^{(r)}_{\text{on}}$, Razo-Mejia et al.
[@Razo-Mejia2018] is used to determine $K_A$, $K_I$, and $k^{(r)}_{\text{on}}$
and Jones et al. [@Jones2014a] is used to determine $r_m$,
$k^{(p)}_{\text{on}}$, and $k^{(p)}_{\text{off}}$.](ch3_fig02){#fig:ch3_fig02
short-caption="Minimal kinetic model of transcriptional regulation for a simple
repression architecture"}
