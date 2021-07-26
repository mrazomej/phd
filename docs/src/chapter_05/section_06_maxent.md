## Maximum entropy approximation of distributions {#sec:ch5_sec06}

(Note: The Python code used for the calculations presented in this section can
be found in the [following
link](https://www.rpgroup.caltech.edu//chann_cap/software/MaxEnt_approx_joint.html)
as an annotated Jupyter notebook)

On the one hand, chemical master equations like the one here represent a hard
mathematical challenge. As presented in Peccoud and Ycart derived a closed-form
solution for the two-state promoter at the mRNA level [@Peccoud1995]. In an
impressive display of mathematical skills, Shahrezaei and Swain were able to
derive an approximate solution for the one- (not considered in this work) and
two-state promoter master equation at the protein level [@Shahrezaei2008].
Nevertheless, both of these solutions do not give instantaneous insights about
the distributions as they involve complicated terms such as confluent
hypergeometric functions.

On the other hand, there has been a great deal of work to generate methods that
can approximate the solution of these discrete state Markovian models [@Ale2013;
@Andreychenko2017; @Frohlich2016; @Schnoerr2017; @Smadbeck2013]. In particular,
for master equations like the one that concerns us here, whose moments can be
easily computed, the moment expansion method provides a simple method to
approximate the full joint distribution of mRNA and protein [@Smadbeck2013].
This section will explain the principles behind this method and show the
implementation for our particular case study.

### The MaxEnt principle

The principle of maximum entropy (MaxEnt) first proposed by E. T. Jaynes in 1957
tackles the question of given limited information what is the least biased
inference one can make about a particular probability distribution
[@Jaynes1957]. In particular, Jaynes used this principle to show the
correspondence between statistical mechanics and information theory,
demonstrating, for example, that the Boltzmann distribution is the probability
distribution that maximizes Shannon's entropy subject to a constraint that the
average energy of the system is fixed.

To illustrate the principle let us focus on a univariate distribution $P_X(x)$.
The $n^{\text{th}}$ moment of the distribution for a discrete set of possible
values of $x$ is given by
$$
    \left\langle{x^n}\right\rangle \equiv \sum_x x^n P_X(x).
    \label{eq:mom_ref}
$$

Now assume that we have knowledge of the first $m$ moments
$\mathbf{\left\langle{x}\right\rangle}_m = (\left\langle{x}\right\rangle,
\left\langle{x^2}\right\rangle, \ldots, \left\langle{x^m}\right\rangle )$. The
question is then how can we use this information to build an estimator $P_H(x
\mid \mathbf{\left\langle{x}\right\rangle}_m)$ of the distribution such that
$$
\lim_{m \rightarrow \infty} 
P_H(x \mid \mathbf{\left\langle{x}\right\rangle}_m) \rightarrow P_X(x),
$$
i.e. that the more moments we add to our approximation, the more the estimator
distribution converges to the real distribution.

The MaxEnt principle tells us that our best guess for this estimator is to build
it based on maximizing the Shannon entropy, constrained by the information we
have about these $m$ moments. Shannon's entropy's maximization guarantees that
we are the least committed possible to information that we do not possess. The
Shannon entropy for a univariate discrete distribution is given by
[@Shannon1948]
$$
H(x) \equiv - \sum_x P_X(x) \log P_X(x).
$$

For an optimization problem subject to constraints, we make use of the method of
the Lagrange multipliers. For this, we define the constraint equation
$\mathcal{L}(x)$ as 
$$
    \mathcal{L}(x) \equiv H(x) - \sum_{i=0}^m
    \left[ \lambda_i \left( \left\langle{x^i}\right\rangle -
    \sum_x x^i P_X(x) \right) \right],
    \label{eq:constraint_eq}
$$
where $\lambda_i$ is the Lagrange multiplier associated with the $i^\text{th}$
moment. The inclusion of the zeroth moment is an additional constraint to
guarantee the normalization of the resulting distribution. Since $P_X(x)$ has a
finite set of discrete values, when taking the derivative of the constraint
equation with respect to $P_X(x)$, we chose a particular value of $X = x$.
Therefore from the sum over all possible $x$ values, only a single term
survives. With this in mind, we take the derivative of the constraint equation,
obtaining 
$$
\frac{d\mathcal{L}}{d P_X(x)} = -\log P_X(x) - 1 -
\sum_{i=0}^m \lambda_i x^i.
$$

Equating this derivative to zero and solving for the distribution (that we now
start calling $P_H(x)$, our MaxEnt estimator) gives
$$
P_H(x) = \exp \left(- 1 - \sum_{i=0}^m \lambda_i x^i \right)
= \frac{1}{\mathcal{Z}}
\exp \left( - \sum_{i=1}^m \lambda_i x^i \right),
\label{eq:maxEnt}
$$
where $\mathcal{Z}$ is the normalization constant that can be obtained by
substituting this solution into the normalization constraint. This results in
$$
\mathcal{Z} \equiv \exp\left( 1 + \lambda_0 \right) =
\sum_x \exp \left( - \sum_{i=1}^m \lambda_i x^i \right).
$$

Eq. $\ref{eq:maxEnt}$ is the general form of the MaxEnt distribution for a
univariate distribution. The computational challenge then consists of finding
numerical values for the Lagrange multipliers $\{ \lambda_i \}$ such that
$P_H(x)$ satisfies our constraints. In other words, the Lagrange multipliers
weigh the contribution of each term in the exponent such that when computing any
of the moments, we recover the value of our constraint. Mathematically what this
means is that $P_H(x)$ must satisfy
$$
\sum_x x^n P_H(x) =
\sum_x \frac{x^n}{\mathcal{Z}}
\exp \left( - \sum_{i=1}^m \lambda_i x^i \right) = 
\left\langle{x^n}\right\rangle.
$$

As an example of applying the MaxEnt principle, let us use a six-face die's
classic problem. If we are only told that after a large number of die rolls, the
mean value of the face is $\left\langle{x}\right\rangle = 4.5$ (note that a fair
die has a mean of $3.5$), what would the least biased guess for the distribution
look like? The MaxEnt principle tells us that our best guess would be of the
form
$$
P_H(x) = \frac{1}{\mathcal{Z}} \exp \left( \lambda x \right).
$$
Using any numerical minimization package, we can easily find the value of the
Lagrange multiplier $\lambda$ that satisfies our constraint. [@Fig:ch5_fig15]
shows two examples of distributions that satisfy the constraint. Panel (A) shows
a distribution consistent with the 4.5 average where both 4 and 5 are equally
likely. Nevertheless, in the information we got about the nature of the die, it
was never stated that some of the faces were forbidden. In that sense, the
distribution is committing to information about the process that we do not
possess. Panel (B), by contrast, shows the MaxEnt distribution that satisfies
this constraint. Since this distribution maximizes Shannon's entropy, it is
guaranteed to be the least biased distribution given the available information.

![**Maximum entropy distribution of six-face die.** (A)biased distribution
consistent with the constraint $\left\langle{x}\right\rangle = 4.5$. (B) MaxEnt
distribution also consistent with the constraint. The Python code
[(`ch5_fig15.py`)](https://github.com/RPGroup-PBoC/chann_cap/blob/master/src/figs/figS15.py)
used to generate this figure can be found on the original paper [GitHub
repository.](https://github.com/RPGroup-PBoC/chann_cap).](ch5_fig15){#fig:ch5_fig15
short-caption="Maximum entropy distribution of six-face die"}

#### The mRNA and protein joint distribution

The MaxEnt principle can easily be extended to multivariate distributions. For
our particular case, we are interested in the mRNA and protein joint
distribution $P(m, p)$. The definition of a moment $\left\langle m^x p^y
\right\rangle$ is a natural extension of Eq. $\ref{eq:mom_ref}$ of the form
$$
\left\langle m^x p^y \right\rangle = \sum_m \sum_p m^x p^y P(m, p).
$$

As a consequence, the MaxEnt joint distribution $P_H(m, p)$ is of the form 
$$
P_H(m, p) = \frac{1}{\mathcal{Z}}
              \exp \left( - \sum_{(x,y)} \lambda_{(x,y)} m^x p^y \right),
\label{eq:maxEnt_joint}
$$
where $\lambda_{(x,y)}$ is the Lagrange multiplier associated with the moment
$\left\langle m^x p^y \right\rangle$, and again $\mathcal{Z}$ is the
normalization constant, given by
$$
\mathcal{Z} = \sum_m \sum_p
              \exp \left( - \sum_{(x, y)} \lambda_{(x, y)} m^x p^y \right).
$$
Note that the sum in the exponent is taken over all available $(x, y)$ pairs
that define the moment constraints for the distribution.

### The Bretthorst rescaling algorithm

The Lagrange multipliers' determination suffers from a numerical underflow and
overflow problem due to the difference in magnitude between the constraints.
This becomes a problem when higher moments are taken into account. The resulting
numerical values for the Lagrange multipliers end up being separated by several
orders of magnitude. For routines such as Newton-Raphson or other minimization
algorithms that can be used to find these Lagrange multipliers, these different
scales become problematic.

To get around this problem, we implemented a variation to the algorithm due to
G. Larry Bretthorst, E.T. Jaynes' last student. With a straightforward argument,
we can show that linearly rescaling the constraints, the Lagrange multipliers,
and the "rules" for computing each of the moments, i.e., each of the individual
products that go into the moment calculation, should converge to the same MaxEnt
distribution. In order to see this let's consider again a univariate
distribution $P_X(x)$ that we are trying to reconstruct given the first two
moments $\left\langle{x}\right\rangle$, and $\left\langle{x^2}\right\rangle$.
The MaxEnt distribution can be written as
$$
P_H(x) = \frac{1}{\mathcal{Z}}
  \exp \left(- \lambda_1 x - \lambda_2 x^2 \right) =
  \frac{1}{\mathcal{Z}}
  \exp \left(- \lambda_1 x \right) \exp \left( - \lambda_2 x^2 \right).
$$
We can always rescale the terms in any way and obtain the same result. Assume
that, for some reason, we want to rescale the quadratic terms by a factor $a$.
We can define a new Lagrange multiplier $\lambda_2' \equiv \frac{\lambda_2}{a}$
that compensates for the rescaling of the terms, obtaining
$$
P_H(x) = \frac{1}{\mathcal{Z}}
  \exp \left(- \lambda_1 x \right) \exp \left( - \lambda_2' ax^2 \right).
$$
Computationally it might be more efficient to find the numerical value of
$\lambda_2'$ rather than $\lambda_2$ maybe because it is of the same order of
magnitude as $\lambda_1$. Then we can always multiply $\lambda_2'$ by $a$ to
obtain back the constraint for our quadratic term. This means that we can always
rescale the MaxEnt problem to make it numerically more stable, then we can
rescale it back to obtain the value of the Lagrange multipliers. The key to the
Bretthorst algorithm lies in selecting what rescaling factor to choose to make
the numerical inference more efficient.

Bretthorst's algorithm goes even further by further transforming the constraints
and the variables to make the constraints orthogonal, making the computation
much more effective. We now explain the algorithm's implementation for our joint
distribution of interest $P(m, p)$.

#### Algorithm implementation

Let the $M \times N$ matrix $\mathbf{A}$ contain all the factors used to compute
the moments that serve as constraints, where each entry is of the form 
$$
A_{ij} = m_i^{x_j} \cdot p_i^{y_j}.
\label{eq:maxent_rules}
$$
In other words, recall that to obtain any moment $\left\langle m^x p^y
\right\rangle$ we compute
$$
\left\langle m^x p^y \right\rangle = \sum_m \sum_p m^x p^y P(m, x).
$$
If we have $M$ possible $(m, p)$ pairs in our truncated sample space (because we
can't include the sample space up to infinity) $\{(m, p)_1, (m, p)_2, \ldots (m,
p)_N \}$, and we have $N$ exponent pairs $(x, y)$ corresponding to the $N$
moments used to constraint the maximum entropy distribution $\{(x, y)_1, (x,
y)_2, \ldots, (x, y)_N \}$, then matrix $\mathbf{A}$ contains all the possible
$M$ by $N$ terms of the form described in Eq. $\ref{eq:maxent_rules}$. Let also
$\mathbf{v}$ be a vector of length $N$ containing all the constraints with each
entry of the form 
$$
v_j = \left\langle{m^{x_j} p^{y_j}}\right\rangle,
$$
i.e. the information that we have about the distribution. That means that the
constraint equation $\mathcal{L}$ to be used for this problem takes the form
$$
\mathcal{L} = -\sum_i P_i \ln P_i + \lambda_0 \left( 1 - \sum_i P_i \right)
  + \sum_{j>0} \lambda_j \left( v_j - \sum_i A_{ij} P_i \right),
$$
where $\lambda_0$ is the Lagrange multiplier associated with the normalization
constraint, and $\lambda_j$ is the Lagrange multiplier associated with the
$j^\text{th}$ constraint. This constraint equation is equivalent to Eq.
$\ref{eq:constraint_eq}$, but now all the details of how to compute the moments
are specified in matrix $\mathbf{A}$.

With this notation in hand, we now proceed to rescale the problem. The first
step consists of rescaling the terms to compute the entries of the matrix
$\mathbf{A}$. As mentioned before, this is the crucial feature of the Bretthorst
algorithm; the particular choice of rescaling factor used in the algorithm
empirically promotes that the rescaled Lagrange multipliers are of the same
order of magnitude. The rescaling takes the form 
$$
A_{ij}' = \frac{A_{ij}}{G_j},
$$
where $G_j$ serves to rescale the moments, providing numerical stability to the
inference problem. Bretthorst proposes an empirical rescaling that satisfies
$$
G_j^2 = \sum_i A_{ij}^2,
$$
or in terms of our particular problem
$$
G_j^2 = \sum_m \sum_p \left( m^{x_j} p^{y_j} \right)^2.
$$
What this indicates is that each pair $m_i^{x_j} p_i^{y_j}$ is normalized by the
square root of the sum of all pairs of the same form squared.

Since we rescale the factors involved in computing the constraints, the
constraints must also be rescaled simply as
$$
v_j' = \left\langle{m^{x_j} p^{y_j}}\right\rangle' = 
\frac{\left\langle{m^{x_j} p^{y_j}}\right\rangle}{G_j}.
$$
The Lagrange multipliers must compensate for this rescaling since the
probability must add up to the same value at the end of the day. Therefore we
rescale the $\lambda_j$ terms as 
$$
\lambda_j' = \lambda_j G_j,
$$
such that any $\lambda_j A_{ij} = \lambda_j' A_{ij}'$. If this empirical value
for the rescaling factor makes the rescaled Lagrange multipliers $\lambda_j'$ be
of the same order of magnitude, this by itself would already improve the
algorithm convergence. Bretthorst proposes another linear transformation to make
the optimization routine even more efficient. For this, we generate orthogonal
constraints that make Newton-Raphson and similar algorithms converge faster. The
transformation is as follows 
$$
A_{ik}'' = \sum_j {e}_{jk} A_{ij}',
$$
for the entires of matrix $\mathbf{A}$, and 
$$
v_k'' = \sum_j {e}_{jk} u_j',
$$
for entires of the constraint vector $\mathbf{v}$, finally 
$$
\lambda_k'' = \sum_j {e}_{jk} \beta_j,
$$
for the Lagrange multipliers. Here ${e}_{jk}$ is the $j^\text{th}$ component of
the $k^\text{th}$ eigenvector of the matrix $\mathbf{E}$ with entries 
$$
{E}_{kj} = \sum_i {A}_{ik}' {A}_{ij}'.
$$
This transformation guarantees that the matrix $\mathbf{A}''$ has the property
$$
\sum_i A_{ij}'' A_{jk}'' = \beta_j \delta_{jk},
$$
where $\beta_j$ is the $j^\text{th}$ eigenvalue of the matrix $\mathbf{E}$ and
$\delta_{jk}$ is the Kronecker delta function. What this means is that, as
desired, the constraints are orthogonal to each other, improving the algorithm
convergence speed.

### Predicting distributions for simple repression constructs

Having explained the theoretical background and the practical difficulties, and
a workaround strategy proposed by Bretthorst, we implemented the inference using
the moments obtained from averaging over the variability along the cell cycle
(See [Sec. 5.4](#sec:ch5_sec05)). [@Fig:ch5_fig16] and [@Fig:ch5_fig17] present these inferences for
both mRNA and protein levels respectively for different values of the
repressor-DNA binding energy and repressor copy numbers per cell. From these
plots, we can easily appreciate that even though the mean of each distribution
changes as the induction level changes, there is a lot of overlap between
distributions. This, as a consequence, means that at the single-cell level,
cells cannot perfectly resolve between different inputs.

![**Maximum entropy mRNA distributions for simple repression constructs.** mRNA
distributions for different biophysical parameters. From left to right, the
repressor-DNA affinity decreases as defined by the three lacI operators O1
($-15.3 \; k_BT$), O2 ($-13.9 \; k_BT$), and O3 ($-9.7 \; k_BT$). From top to
bottom, the mean repressor copy number per cell increases. The curves on each
plot represent different IPTG concentrations. Each distribution was fitted using
the first three moments of the mRNA distribution. The Python code
[(`ch5_fig16.py`)](https://github.com/RPGroup-PBoC/chann_cap/blob/master/src/figs/figS16.py)
used to generate this figure can be found on the original paper [GitHub
repository.](https://github.com/RPGroup-PBoC/chann_cap).](ch5_fig16){#fig:ch5_fig16
short-caption="Maximum entropy mRNA distributions for simple repression
constructs"}

![**Maximum entropy protein distributions for simple repression constructs.**
Protein distributions for different biophysical parameters. From left to right,
the repressor-DNA affinity decreases as defined by the three lacI operators O1
($-15.3 \; k_BT$), O2 ($-13.9 \; k_BT$), and O3 ($-9.7 \; k_BT$). From top to
bottom, the mean repressor copy number per cell increases. The curves on each
plot represent different IPTG concentrations. Each distribution was fitted using
the first six moments of the protein distribution. The Python code
[(`ch5_fig17.py`)](https://github.com/RPGroup-PBoC/chann_cap/blob/master/src/figs/figS17.py)
used to generate this figure can be found on the original paper [GitHub
repository.](https://github.com/RPGroup-PBoC/chann_cap).](ch5_fig17){#fig:ch5_fig17
short-caption="Maximum entropy protein distributions for simple repression
constructs"}

### Comparison with experimental data

Now that we have reconstructed an approximation of the probability distribution
$P(m, p)$ we can compare this with our experimental measurements. But just as
detailed in the single-cell microscopy, measurements are given in arbitrary
units of fluorescence. Therefore we cannot directly compare our predicted
protein distributions with these values. To get around this issue, we use the
fact that the fold-change in gene expression that we defined as the ratio of the
gene expression level in the presence of the repressor and the expression level
of a knockout strain is a non-dimensional quantity. Therefore we normalize all
of our single-cell measurements by the mean fluorescence value of the $\Delta
lacI$ strain with the proper background fluorescence subtracted as explained in
the noise measurements. In the case of the theoretical predictions of the
protein distribution, we also normalize each protein value by the predicted mean
protein level $\left\langle p \right\rangle$, having now non-dimensional scales
that can be directly compared. [@Fig:ch5_fig18] shows the experimental (color
curves) and theoretical (dark dashed line) cumulative distribution functions for
the three $\Delta lacI$ strains. As in [@Fig:ch5_fig12], we do not expect
differences between the operators, but we explicitly plot them separately to
ensure that this is the case. We can see right away that as we would expect,
given the limitations of the model to predict the noise and skewness of the
distribution accurately, the model doesn't accurately predict the data. Our
model predicts a narrower distribution compared to what we measured with
single-cell microscopy.

![**Experiment vs. theory comparison for $\Delta lacI$ strain.** Example
fold-change empirical cumulative distribution functions (ECDF) for strains with
no repressors and different operators. The color curves represent single-cell
microscopy measurements while the dashed black lines represent the theoretical
distributions as reconstructed by the maximum entropy principle. The theoretical
distributions were fitted using the first six moments of the protein
distribution. The Python code
[(`ch5_fig18.py`)](https://github.com/RPGroup-PBoC/chann_cap/blob/master/src/figs/figS18.py)
used to generate this figure can be found on the original paper [GitHub
repository.](https://github.com/RPGroup-PBoC/chann_cap).](ch5_fig18){#fig:ch5_fig18
short-caption="Experiment vs. theory comparison for $\Delta lacI$ strain"}

The same narrower prediction applies to the regulated promoters.
[@Fig:ch5_fig19], shows the theory-experiment comparison of the cumulative
distribution functions for different repressor binding sites (different
figures), repressor copy numbers (rows), and inducer concentrations (columns).
In general, the predictions are systematically narrower compared to the actual
experimental data.

![**Experiment vs. theory comparison for regulated promoters.** Example
fold-change empirical cumulative distribution functions (ECDF) for regulated
strains with the three operators (different colors) as a function of repressor
copy numbers (rows) and inducer concentrations (columns). The color curves
represent single-cell microscopy measurements while the dashed black lines
represent the theoretical distributions as reconstructed by the maximum entropy
principle. The theoretical distributions were fitted using the first six moments
of the protein distribution. The Python code
[(`ch5_fig19.py`)](https://github.com/RPGroup-PBoC/chann_cap/blob/master/src/figs/figS19.py)
used to generate this figure can be found on the original paper [GitHub
repository.](https://github.com/RPGroup-PBoC/chann_cap).](ch5_fig19){#fig:ch5_fig19
short-caption="Experiment vs. theory comparison for regulated promoters"}
