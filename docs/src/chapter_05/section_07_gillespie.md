## Gillespie simulation of the master equation 

(Note: The Python code used for the calculations presented in this section can
be found in the [following
link](https://www.rpgroup.caltech.edu//chann_cap/software/gillespie_simulation.html)
as an annotated Jupyter notebook)

So far, we have generated a way to compute an approximated form of the joint
distribution of protein and mRNA $P(m, p)$ as a function of the moments of the
distribution $\left\langle m^x p^y \right\rangle$. This is a non-conventional
form to work with the resulting distribution of the master equation. A more
conventional approach to work with master equations whose closed-form solutions
are not known or not computable is to use stochastic simulations, commonly known
as Gillespie simulations. To benchmark our approach's performance based on
distribution moments and maximum entropy, we implemented the Gillespie
algorithm. Our implementation, as detailed in the corresponding Jupyter
notebook, makes use of just-in-time compilation as implemented with the Python
package [numba](http://numba.pydata.org).

### mRNA distribution with Gillespie simulations

To confirm that the Gillespie simulation's implementation was correct, we
perform the simulation at the mRNA level, for which the closed-form solution of
the steady-state distribution is known as detailed in [Sec.
5.2](#sec:ch5_sec03). [@Fig:ch5_fig20] shows example trajectories of mRNA
counts. Each of these trajectories was computed over several cell cycles, where
the cell division was implemented, generating a binomially distributed random
variable that depended on the last mRNA count before the division event.

![**Stochastic trajectories of mRNA counts.** 100 stochastic trajectories
generated with the Gillespie algorithm for mRNA counts over time for a two-state
unregulated promoter. Cells spend a fraction of the cell cycle with a single
copy of the promoter (light brown) and the rest of the cell cycle with two
copies (light yellow). When trajectories reach a new cell cycle, the mRNA counts
undergo binomial partitioning to simulate the cell division. The Python code
[(`ch5_fig20.py`)](https://github.com/RPGroup-PBoC/chann_cap/blob/master/src/figs/figS20.py)
used to generate this figure can be found on the original paper [GitHub
repository.](https://github.com/RPGroup-PBoC/chann_cap).](ch5_fig20){#fig:ch5_fig20
short-caption="Stochastic trajectories of mRNA counts"}

To check the implementation of our stochastic algorithm, we generated several of
these stochastic trajectories to reconstruct the mRNA steady-state distribution.
These reconstructed distributions for a single- and double-copy of the promoter
can be compared with Eq. $\ref{eq:two_state_mRNA}$--the steady-state
distribution for the two-state promoter. [@Fig:ch5_fig21] shows the excellent
agreement between the stochastic simulation and the analytical result,
confirming that our implementation of the Gillespie simulation is correct.

![**Comparison of analytical and simulated mRNA distribution.** Solid lines show
the steady-state mRNA distributions for one copy (light blue) and two copies of
the promoter (dark blue) as defined by Eq. $\ref{eq:two_state_mRNA}$. Shaded
regions represent the corresponding distribution obtained using 2500 stochastic
mRNA trajectories and taking the last cell cycle to approximate the
distribution. The Python code
[(`ch5_fig21.py`)](https://github.com/RPGroup-PBoC/chann_cap/blob/master/src/figs/figS21.py)
used to generate this figure can be found on the original paper [GitHub
repository.](https://github.com/RPGroup-PBoC/chann_cap).](ch5_fig21){#fig:ch5_fig21
short-caption="Comparison of analytical and simulated mRNA distribution"}

### Protein distribution with Gillespie simulations

Having confirmed that our implementation of the Gillespie algorithm that
includes the binomial partitioning of molecules reproduces analytical results,
we extended the implementation to include protein counts. [@Fig:ch5_fig22] shows
representative trajectories for both mRNA and protein counts over several cell
cycles. Especially for the protein, we can see that it takes several cell cycles
for counts to converge to the dynamical steady-state observed with the
deterministic moment equations. Once this steady-state is reached, the ensemble
of trajectories between cell cycles looks very similar.

![**Stochastic trajectories of mRNA and protein counts.** 2500 protein counts
over time for a two-state unregulated promoter. Cells spend a fraction of the
cell cycle with a single copy of the promoter (light brown) and the rest of the
cell cycle with two copies (light yellow). When trajectories reach a new cell
cycle, the molecule counts undergo binomial partitioning to simulate the cell
division. The Python code
[(`ch5_fig22.py`)](https://github.com/RPGroup-PBoC/chann_cap/blob/master/src/figs/figS22.py)
used to generate this figure can be found on the original paper [GitHub
repository.](https://github.com/RPGroup-PBoC/chann_cap).](ch5_fig22){#fig:ch5_fig22
short-caption="Stochastic trajectories of mRNA and protein counts"}

From these trajectories, we can compute the steady-state protein distribution,
taking into account the cell-age distribution, as detailed in [Sec.
5.5](#sec:ch5_sec06). [@Fig:ch5_fig23] shows the comparison between this
distribution and the one generated using the maximum entropy algorithm. Although
the notorious differences between the distributions, the Gillespie simulation
and the maximum entropy results are indistinguishable in terms of the mean,
variance, and skewness of the distribution. We remind the reader that the
maximum entropy approximates the distribution that gets better the more moments
we add. We, therefore, claim that the approximation works sufficiently well for
our purpose. The enormous advantage of the maximum entropy approach comes from
the computation time. For the number of distributions needed for our
calculations, the Gillespie algorithm proved to be a very inefficient method
given the ample sample space. Our maximum entropy approach reduces the
computation time by several orders of magnitude, allowing us to explore
different regulatory models' parameters extensively.

![**Comparison of protein distributions.** Comparison of the protein
distribution generated with Gillespie stochastic simulations (blue curve) and
the maximum entropy approach (orange curve). The upper panel shows the
probability mass function. The lower panel compares the cumulative distribution
functions. The Python code
[(`ch5_fig23.py`)](https://github.com/RPGroup-PBoC/chann_cap/blob/master/src/figs/figS23.py)
used to generate this figure can be found on the original paper [GitHub
repository.](https://github.com/RPGroup-PBoC/chann_cap).](ch5_fig23){#fig:ch5_fig23
short-caption="Comparison of protein distributions"}
