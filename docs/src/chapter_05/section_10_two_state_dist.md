## Derivation of the steady-state mRNA distribution

In this section we will derive the two-state promoter mRNA distribution we quote
in XXX. For this method we will make use of the so-called generating functions.
Generating functions are mathematical objects on which we can encode a series of
infinite numbers as coefficients of a power series. The power of generating
functions comes from the fact that we can convert an infinite-dimensional system
of coupled ordinary differential equations--in our case the system of
differential equations defining all probabilities $P(m, t)$ for $m \in
\mathbb{Z}$--into a single partial differential equation that we can then solve
to extract back the probability distributions.

To motivate the use of generating functions we will begin with the simplest
case: the one-state Poisson promoter.

### One-state Poisson promoter

We begin by defining the reaction scheme that defines the one-state promoter.
[@Fig:ch5_fig35] shows the schematic representation of the Poisson promoter as a
simple cartoon (part (A)) and as the Markov chain that defines the state space
of the system (part (B)).

![**One-state Poisson promoter.** (A) Schematic of the kinetics of the one
state-promoter. mRNA is produced and degrade stochastically with a rate $r_m$
and $\gamma_m$, respectively. (B) Representation of the Markov-chain for the
state space that the promoter can be in. The distribution $P(m, t)$ represents
the probability of having certain discrete number of mRNA $m$ at time $t$. The
transition between states depends on the previously mentioned
rates.](ch5_fig35){#fig:ch5_fig35 short-caption="One-state Poisson promoter"}

The dynamics of the probability distribution $P(m, t)$ are governed by the
chemical master equation
$$
\frac{d P(m, t)}{dt} = 
\overbrace{r_m P(m-1, t)}^{m-1 \rightarrow m}
- \overbrace{r_m P(m, t)}^{m \rightarrow m+1}
+ \overbrace{\gamma_m (m+1) P(m+1, t)}^{m+1 \rightarrow m}
- \overbrace{\gamma_m m P(m, t)}^{m \rightarrow m-1}.
\label{eq:eq_test}
$$

This is $\ref{eq:eq_test}$
