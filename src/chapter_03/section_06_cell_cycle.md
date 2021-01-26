### Accounting for cell-cycle dependent variability in gene dosage

As cells progress through the cell cycle, the genome has to be replicated to
guarantee that each daughter cell receives a copy of the genetic material. As
replication of the genome can take longer than the total cell cycle, this
implies that cells spend part of the cell cycle with multiple copies of each
gene depending on the cellular growth rate and the relative position of the gene
with respect to the replication origin [@Bremer1996]. Genes closer to the
replication origin spend a larger fraction of the cell cycle with multiple
copies compared to genes closer to the replication termination site
[@Bremer1996]. (A) depicts a schematic of this process where the replication
origin (*oriC*) and the relevant locus for our experimental measurements
(*galK*) are highlighted.

Since this change in gene copy number has been shown to have an effect on
cell-to-cell variability in gene expression [@Jones2014a; @Peterson2015], we now
extend our minimal model to account for these changes in gene copy number during
the cell cycle. We reason that the only difference between the single-copy state
and the two-copy state of the promoter is a doubling of the mRNA production rate
$r_m$. In particular, the promoter activation and inactivation rates
$k^{(p)}_{\text{on}}$ and $k^{(p)}_{\text{off}}$ and the mRNA production rate
$r_m$ inferred in assume that cells spend a fraction $f$ of the cell cycle with
one copy of the promoter (mRNA production rate $r_m$) and a fraction $(1-f)$ of
the cell cycle with two copies of the promoter (mRNA production rate $2 r_m$).
This inference was performed considering that at each cell state the mRNA level
immediately reaches the steady state value for the corresponding mRNA production
rate. This assumption is justified since the timescale to reach this steady
state depends only on the degradation rate $\gamma _m$, which for the mRNA is
much shorter ($\approx 3$ min) than the length of the cell cycle ($\approx$ 60
min for our experimental conditions) [@Dong1995]. shows that a model accounting
for this gene copy number variability is able to capture data from single
molecule mRNA counts of an unregulated (constitutively expressed) promoter.

Given that the protein degradation rate $\gamma _p$ in our model is set by the
cell division time, we do not expect that the protein count will reach the
corresponding steady state value for each stage in the cell cycle. In other
words, cells do not spend long enough with two copies of the promoter for the
protein level to reach the steady state value corresponding to a transcription
rate of $2 r_m$. We therefore use the dynamical equations developed in to
numerically integrate the time trajectory of the moments of the distribution
with the corresponding parameters for each phase of the cell cycle. (B) shows an
example corresponding to the mean mRNA level (upper panel) and the mean protein
level (lower panel) for the case of the unregulated promoter. Given that we
inferred the promoter rate parameters considering that mRNA reaches steady state
in each stage, we see that the numerical integration of the equations is
consistent with the assumption of having the mRNA reach a stable value in each
stage (See (B) upper panel). On the other hand, the mean protein level does not
reach a steady state at either of the cellular stages. Nevertheless it is
notable that after several cell cycles the trajectory from cycle to cycle
follows a repetitive pattern (See (B) lower panel). Previously we have
experimentally observed this repetitive pattern by tracking the expression level
over time with video microscopy as observed in Fig. 18 of [@Phillips2019].

To test the effects of including this gene copy number variability in our model
we now compare the predictions of the model with experimental data. As detailed
in the Methods section, we obtained single-cell fluorescence values of different
*E. coli* strains carrying a YFP gene under the control of the LacI repressor.
Each strain was exposed to twelve different input inducer (IPTG) concentrations
for $\approx 8$ generations for cells to adapt to the media. The strains imaged
spanned three orders of magnitude in repressor copy number and three distinct
repressor-DNA affinities. Since growth was asynchronous, we reason that cells
were randomly sampled at all stages of the cell cycle. Therefore, when computing
statistics from the data such as the mean fluorescence value, in reality we are
averaging over the cell cycle. In other words, as depicted in (B), quantities
such as the mean protein copy number change over time, i.e. $\langle p \rangle
\equiv \langle p(t) \rangle$. This means that computing the mean of a population
of unsynchronized cells is equivalent to averaging this time-dependent mean
protein copy number over the span of the cell cycle. Mathematically this is
expressed as
$$
\langle p \rangle_c = \int_{t_o}^{t_d} \langle p(t) \rangle P(t) dt,
$${#eq:ch3_eq13}
where $\langle p(t) \rangle$ represents the first moment of the protein
distribution as computed from , $\langle p\rangle_c$ represents the average
protein copy number over a cell cycle, $t_o$ represents the start of the cell
cycle, $t_d$ represents the time of cell division, and $P(t)$ represents the
probability of any cell being at time $t \in [t_o, t_d]$ of their cell cycle. We
do not consider cells uniformly distributed along the cell cycle since it is
known that cells age is exponentially distributed, having more younger than
older cells at any point in time [@Powell1956] (See for further details). All
computations hereafter are therefore done by applying an average like that in
for the span of a cell cycle. We remind the reader that these time averages are
done under a fixed environmental state. It is the trajectory of cells over cell
cycles under a constant environment that we need to account for. It is through
this averaging over the span of a cell cycle that we turn a periodic process as
the one shown in [@Fig:ch3_fig03](B) into a stationary process that we can
compare with experimental data and, as we will see later, use to reconstruct the
steady state gene expression distribution.

[@Fig:ch3_fig03](C) compares zero-parameter fit predictions (lines) with
experimentally determined quantities (points). The upper row shows the
non-dimensional quantity known as the fold-change in gene expression
[@Garcia2011c]. This fold-change is defined as the relative mean gene expression
level with respect to an unregulated promoter. For protein this is 
$$
\text{fold-change} = 
\frac{\langle p(R > 0) \rangle_c}{\langle p(R = 0) \rangle_c},
$${#eq:ch2_eq14}
where $\langle p(R > 0)i \rangle_c$ represents the mean protein count for cells
with non-zero repressor copy number count $R$ over the entire cell cycle, and
$\langle p(R = 0) \rangle_c$ represents the equivalent for a strain with no
repressors present. The experimental points were determined from the YFP
fluorescent intensities of cells with varying repressor copy number and a
$\Delta lacI$ strain with no repressor gene present (See Methods for further
details). The fold-change in gene expression has previously served as a metric
to test the validity of equilibrium-based models [@Phillips2015a]. We note that
the curves shown in the upper panel of (C) are consistent with the predictions
from equilibrium models [@Razo-Mejia2018] despite being generated from a clearly
non-equilibrium process as shown in (B). The kinetic model from (A) goes beyond
the equilibrium picture to generate predictions for moments of the distribution
other than the mean mRNA or mean protein count. To test this extended predictive
power the lower row of (C) shows the noise in gene expression defined as the
standard deviation over the mean protein count, accounting for the changes in
gene dosage during the cell cycle. Although our model systematically
underestimates the noise in gene expression, the zero-parameter fits capture the
scaling of this noise. Possible origins of this systematic discrepancy could be
the intrinsic cell-to-cell variability of rate parameters given the variability
in the molecular components of the central dogma machinery [@Jones2014a], or
noise generated by irreversible non-equilibrium reactions not explicitly taken
into account in our minimal model [@Grah2020]. The large errors for the highly
repressed strains (lower left panel in (C)) are a result of having a small
number in the denominator - mean fluorescence level - when computing the noise.
Although the model is still highly informative about the physical nature of how
cells regulate their gene expression, the lack of exact numerical agreement
between theory and data opens an opportunity to gain new insights into the
biophysical origin of cell-to-cell variability. In we explore empirical ways to
account for this systematic deviation. We direct the reader to where equivalent
predictions are done ignoring the changes in gene dosage due to the replication
of the genome.

![**Accounting for gene copy number variability during the cell cycle.** (A)
Schematic of a replicating bacterial genome. As cells progress through the cell
cycle the genome is replicated, duplicating gene copies for a fraction of the
cell cycle before the cell divides. *oriC* indicates the replication origin, and
*galK* indicates the locus at which the YFP reporter construct was integrated.
(B) mean (solid line) $\pm$ standard deviation (shaded region) for the mRNA
(upper panel) and protein (lower panel) dynamics. Cells spend a fraction of the
cell cycle with a single copy of the promoter (light brown) and the rest of the
cell cycle with two copies (light yellow). Black arrows indicate time of cell
division. (C) Zero parameter-fit predictions (lines) and experimental data
(circles) of the gene expression fold-change (upper row) and noise (lower row)
for repressor binding sites with different affinities (different columns) and
different repressor copy numbers per cell (different lines on each panel). Error
bars in data represent the 95% confidence interval on the quantities as computed
from 10,000 bootstrap estimates generated from $> 500$ single-cell fluorescence
measurements. In the theory curves, dotted lines indicate plot in linear scale
to include zero, while solid lines indicate logarithmic scale. For visual
clarity, data points in the noise panel with exceptionally large values coming
from highly repressed strains are plotted on a separate
panel.](ch3_fig03){#fig:ch3_fig03 short-caption="Accounting for gene copy number
variability during the cell cycle."}