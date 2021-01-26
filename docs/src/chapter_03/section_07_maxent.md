### Maximum Entropy approximation

Having numerically computed the moments of the mRNA and protein distributions as
cells progress through the cell cycle, we now proceed to make an approximate
reconstruction of the full distributions given this limited information. As
hinted in the maximum entropy principle, first proposed by E.T. Jaynes in 1957
[@Jaynes1957], approximates the entire distribution by maximizing the Shannon
entropy subject to constraints given by the values of the moments of the
distribution [@Jaynes1957]. This procedure leads to a probability distribution
of the form (See for full derivation)
$$
P(m, p) = \frac{1}{\mathcal{Z}} \exp
\left( - \sum_{(x,y)} \lambda_{(x,y)} m^x p^y \right), 
$${#eq:ch3_eq15}
where $\lambda_{(x,y)}$ is the Lagrange multiplier associated with the
constraint set by the moment $\langle m^x p^y \rangle$, and $\mathcal{Z}$ is a
normalization constant. The more moments $\langle m^x p^y \rangle$ included as
constraints, the more accurate the approximation resulting from becomes.

The computational challenge then becomes an optimization routine in which the
values for the Lagrange multipliers $\lambda_{(x,y)}$ that are consistent with
the constraints set by the moment values $\langle m^x p^y \rangle$ need to be
found. This is computationally more efficient than sampling directly out of the
master equation with a stochastic algorithm (see for further comparison between
maximum entropy estimates and the Gillespie algorithm). details our
implementation of a robust algorithm to find the values of the Lagrange
multipliers. (A) shows example predicted protein distributions reconstructed
using the first six moments of the protein distribution for a suite of different
biophysical parameters and environmental inducer concentrations. As
repressor-DNA binding affinity (columns in (A)) and repressor copy number (rows
in (A)) are varied, the responses to different signals, i.e. inducer
concentrations, overlap to varying degrees. For example, the upper right corner
frame with a weak binding site ($\Delta\varepsilon_r = -9.7 \; k_BT$) and a low
repressor copy number (22 repressors per cell) have virtually identical
distributions regardless of the input inducer concentration. This means that
cells with this set of parameters cannot resolve any difference in the
concentration of the signal. As the number of repressors is increased, the
degree of overlap between distributions decreases, allowing cells to better
resolve the value of the signal input. On the opposite extreme the lower left
panel shows a strong binding site ($\Delta\varepsilon_r = -15.3 \; k_BT$) and a
high repressor copy number (1740 repressors per cell). This parameter
combination shows overlap between distributions since the high degree of
repression centers all distributions towards lower copy numbers, again giving
little ability for the cells to resolve the inputs. In (B) and we show the
comparison of these predicted cumulative distributions with the experimental
single-cell fluorescence distributions. Given the systematic deviation of our
predictions for the protein copy number noise highlighted in (C), the
theoretical distributions (dashed lines) underestimate the width of the
experimental data. We again direct the reader to for an exploration of empirical
changes to the moments that improve the agreement of the predictions. In the
following section we formalize the notion of how well cells can resolve
different inputs from an information theoretic perspective via the channel
capacity.

![**Maximum entropy protein distributions for varying physical parameters.** (A)
Predicted protein distributions under different inducer (IPTG) concentrations
for different combinations of repressor-DNA affinities (columns) and repressor
copy numbers (rows). The first six moments of the protein distribution used to
constrain the maximum entropy approximation were computed by integrating as
cells progressed through the cell cycle as described in . (B) Theory-experiment
comparison of predicted fold-change empirical cumulative distribution functions
(ECDF). Each panel shows two example concentrations of inducer (colored curves)
with their corresponding theoretical predictions (dashed lines). Distributions
were normalized to the mean expression value of the unregulated strain in order
to compare theoretical predictions in discrete protein counts with experimental
fluorescent measurements in arbitrary units.](ch3_fig04){#fig:ch3_fig04
short-caption="Maximum entropy protein distributions for varying physical
parameters."}