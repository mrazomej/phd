## Computing moments from the master equation

In this section we will compute the moment equations for the distribution $P(m,
p)$. Without lost of generality here we will focus on the three-state regulated
promoter. The computation of the moments for the two-state promoter follows the
exact same procedure, changing only the definition of the matrices in the master
equation.

### Computing moments of a distribution

(Note: The Python code used for the calculations presented in this section can
be found in the [following
link](https://www.rpgroup.caltech.edu//chann_cap/software/moment_dynamics_system.html)
as an annotated Jupyter notebook)

To compute any moment of our chemical master equation () let us define a
vector
\begin{equation}
\left\langle \mathbf{m^x p^y} \right\rangle \equiv
(\left\langle m^x p^y \right\rangle_A, 
\left\langle m^x p^y \right\rangle_I, \left\langle m^x p^y \right\rangle_R)^T,
\end{equation}
where $\left\langle m^x p^y \right\rangle_S$ is the expected value of $m^x p^y$
in state $S \in \{A, I, R\}$ with $x, y \in \mathbb{N}$. In other words, just as
we defined the vector $\mathbf{P}(m, p)$, here we define a vector to collect the
expected value of each of the promoter states. By definition, these moments
$\left\langle m^x p^y \right\rangle_S$ are computed as
\begin{equation}
\left\langle m^x p^y \right\rangle_S \equiv 
\sum_{m=0}^\infty \sum_{p=0}^\infty m^x p^y P_S(m, p).
\end{equation}
To simplify the notation, let $\sum_x \equiv \sum_{x=0}^\infty$. Since we are
working with a system of three ODEs, one for each state, let us define the
following operation:
\begin{equation}
\left\langle \mathbf{m^x p^y} \right\rangle =
\sum_m \sum_p m^x p^y \mathbf{P}(m, p) \equiv
  \begin{bmatrix}
    \sum_m \sum_p m^x p^y P_A(m, p)\\
    \sum_m \sum_p m^x p^y P_I(m, p)\\
    \sum_m \sum_p m^x p^y P_R(m, p)\\
  \end{bmatrix}.
\end{equation}

With this in hand we can then apply this sum over $m$ and $p$ to . For the
left-hand side we have
\begin{equation}
\sum_m \sum_p m^x p^y \frac{d \mathbf{P}(m, p)}{dt} = 
\frac{d}{dt}\left[ \sum_m \sum_p m^x p^y \mathbf{P}(m, p) \right],
\end{equation}
where we made use of the linearity property of the derivative to switch the
order between the sum and the derivative. Notice that the right-hand side of
contains the definition of a moment from XXX. That means that we can rewrite it
as
\begin{equation}
\frac{d}{dt}\left[ \sum_m \sum_p m^x p^y \mathbf{P}(m, p) \right] = 
\frac{d \mathbf{\left\langle m^x p^y \right\rangle}}{dt}.
\end{equation}

Distributing the sum on the right-hand side of gives 
\begin{equation}
\begin{split}
    \frac{d \mathbf{\left\langle m^x p^y \right\rangle}}{dt} &=
    \mathbf{K} \sum_m \sum_p m^x p^y \mathbf{P}(m, p)\\
    &- \mathbf{R}_m \sum_m \sum_p m^x p^y \mathbf{P}(m, p) + 
    \mathbf{R}_m \sum_m \sum_p m^x p^y \mathbf{P}(m-1, p)\\
    &- \mathbf{\Gamma}_m \sum_m \sum_p (m) m^x p^y \mathbf{P}(m, p) + 
    \mathbf{\Gamma}_m \sum_m \sum_p (m + 1) m^x p^y \mathbf{P}(m + 1, p)\\
    &- \mathbf{R}_p \sum_m \sum_p (m) m^x p^y \mathbf{P}(m, p) + 
    \mathbf{R}_p \sum_m \sum_p (m) m^x p^y \mathbf{P}(m, p - 1)\\
    &- \mathbf{\Gamma}_p \sum_m \sum_p (p) m^x p^y \mathbf{P}(m, p) + 
    \mathbf{\Gamma}_p \sum_m \sum_p (p + 1) m^x p^y \mathbf{P}(m, p + 1).
  \end{split}
\end{equation}

Let's look at each term on the right-hand side individually. For the terms in
involving $\mathbf{P}(m, p)$ we can again use to rewrite them in a more compact
form. This means that we can rewrite the state transition term as 
\begin{equation}
\mathbf{K} \sum_m \sum_p m^x p^y \mathbf{P}(m, p) = 
\mathbf{K} \mathbf{\left\langle m^x p^y \right\rangle}.
\end{equation}
The mRNA production term involving $\mathbf{P}(m, p)$ can be rewritten as
\begin{equation}
\mathbf{R}_m \sum_m \sum_p m^x p^y \mathbf{P}(m, p) = 
\mathbf{R}_m \mathbf{\left\langle m^x p^y \right\rangle}.
\end{equation}
In the same way the mRNA degradation term gives
\begin{equation}
\mathbf{\Gamma}_m \sum_m \sum_p (m) m^x p^y \mathbf{P}(m, p) = 
\mathbf{\Gamma}_m \mathbf{\left\langle{m^{(x + 1)} p^y}\right\rangle}.
\end{equation}
For the protein production and degradation terms involving $\mathbf{P}(m, p)$ we
have 
\begin{equation}
\mathbf{R}_p \sum_m \sum_p (m) m^x p^y \mathbf{P}(m, p) = 
\mathbf{R}_p \mathbf{\left\langle{m^{(x + 1)} p^y}\right\rangle},
\end{equation}
and 
\begin{equation}
\mathbf{\Gamma}_p \sum_m \sum_p (p) m^x p^y \mathbf{P}(m, p) = 
\mathbf{\Gamma}_p \mathbf{\left\langle{m^x p^{(y + 1)}}\right\rangle},
\end{equation}
respectively.

For the sums terms in involving $\mathbf{P}(m \pm 1, p)$ or $\mathbf{P}(m, p \pm
1)$ we can reindex the sum to work around this mismatch. To be more specific
let's again look at each term case by case. For the mRNA production term
involving $\mathbf{P}(m-1, p)$ we define $m' \equiv m - 1$. Using this we write
\begin{equation}
\mathbf{R}_m \sum_m \sum_p m^x p^y \mathbf{P}(m-1, p) =
\mathbf{R}_m \sum_{m' = -1}^\infty \sum_p (m' + 1)^x p^y \mathbf{P}(m', p).
\end{equation}
Since having negative numbers of mRNA or protein doesn't make physical sense we
have that $\mathbf{P}(-1, p) = 0$. Therefore we can rewrite the sum starting
from 0 rather than from -1, obtaining
\begin{equation}
\mathbf{R}_m \sum_{m' = -1}^\infty \sum_p (m' + 1)^x p^y \mathbf{P}(m', p) =
\mathbf{R}_m \sum_{m'=0}^\infty \sum_p (m' + 1)^x p^y \mathbf{P}(m', p).
\end{equation}
Recall that our distribution $\mathbf{P}(m, p)$ takes $m$ and $p$ as numerical
inputs and returns a probability associated with such a molecule count.
Nevertheless, $m$ and $p$ themselves are dimensionless quantities that serve as
indices of how many molecules are in the cell. The distribution is the same
whether the variable is called $m$ or $m'$; for a specific number, let's say $m
= 5$, or $m' = 5$, $\mathbf{P}(5, p)$ will return the same result. This means
that the variable name is arbitrary, and the right-hand side of can be written
as
\begin{equation}
\mathbf{R}_m \sum_{m'=0}^\infty \sum_p (m' + 1)^x p^y \mathbf{P}(m', p) =
\mathbf{R}_m \mathbf{\left\langle{(m+1)^x p^y}\right\rangle},
\end{equation}
since the left-hand side corresponds to the definition of a moment.

For the mRNA degradation term involving $\mathbf{P}(m + 1, p)$ we follow a
similar procedure in which we define $m' = m + 1$ to obtain
\begin{equation}
\mathbf{\Gamma}_m \sum_m \sum_p (m + 1) m^x p^y \mathbf{P}(m + 1, p) =
\mathbf{\Gamma}_m \sum_{m' = 1}^\infty \sum_p m' 
(m' - 1)^x p^y \mathbf{P}(m', p).
\end{equation}
In this case since the term on the right-hand side of the equation is multiplied
by $m'$, starting the sum over $m'$ from 0 rather than from 1 will not affect
the result since this factor will not contribute to the total sum. Nevertheless
this is useful since our definition of a moment from requires the sum to start
at zero. This means that we can rewrite this term as
\begin{equation}
\mathbf{\Gamma}_m \sum_{m' = 1}^\infty m' \sum_p 
(m' - 1)^x p^y \mathbf{P}(m', p) =
\mathbf{\Gamma}_m \sum_{m' = 0}^\infty m' \sum_p 
(m' - 1)^x p^y \mathbf{P}(m', p).
\end{equation}
Here again we can change the arbitrary label $m'$ back to $m$ obtaining
\begin{equation}
\mathbf{\Gamma}_m \sum_{m' = 0}^\infty m' \sum_p 
(m' - 1)^x p^y \mathbf{P}(m', p) =
\mathbf{\Gamma}_m \mathbf{\left\langle{m (m - 1)^x p^y}\right\rangle}.
\end{equation}

The protein production term involving $\mathbf{P}(m, p - 1)$ can be reindexed by
defining $p' \equiv p - 1$. This gives
\begin{equation}
\mathbf{R}_p \sum_m \sum_p (m) m^x p^y \mathbf{P}(m, p - 1) =
\mathbf{R}_p \sum_m \sum_{p'=-1}^\infty 
m^{(x + 1)} (p + 1)^y \mathbf{P}(m, p').
\end{equation}
We again use the fact that negative molecule copy numbers are assigned with
probability zero to begin the sum from 0 rather than -1 and the arbitrary nature
of the label $p'$ to write
\begin{equation}
\mathbf{R}_p \sum_m \sum_{p'=0}^\infty m^{(x + 1)} (p + 1)^y \mathbf{P}(m, p') =
\mathbf{R}_p \mathbf{\left\langle{m^{(x + 1)} (p + 1)^y}\right\rangle}.
\end{equation}
Finally, we take care of the protein degradation term involving $\mathbf{P}(m, p
+ 1)$. As before, we define $p' = p + 1$ and substitute this to obtain
\begin{equation}
\mathbf{\Gamma}_p \sum_m \sum_p (p + 1) m^x p^y \mathbf{P}(m, p + 1) =
\mathbf{\Gamma}_p \sum_m \sum_{p'=1}^\infty 
(p') m^x (p' - 1)^y \mathbf{P}(m, p').
\end{equation}
Just as with the mRNA degradation term, having a term $p'$ inside the sum allows
us to start the sum over $p'$ from 0 rather than 1. We can therefore write
\begin{equation}
\mathbf{\Gamma}_p \sum_m \sum_{p'=0}^\infty 
(p') m^x (p' - 1)^y \mathbf{P}(m, p') =
\mathbf{\Gamma}_p \mathbf{\left\langle{m^x p (p - 1)^y}\right\rangle}.
\end{equation}

Putting all these terms together we can write the general moment ODE. This is of
the form 
\begin{equation}
\begin{split}
    \frac{d\mathbf{\left\langle m^x p^y \right\rangle}}{dt} &=
    \mathbf{K} \mathbf{\left\langle m^x p^y \right\rangle}
    \text{  (promoter state transition)}\\
    &- \mathbf{R}_m \mathbf{\left\langle m^x p^y \right\rangle} +
    \mathbf{R}_m \mathbf{\left\langle{(m+1)^x p^y}\right\rangle}
    \text{  (mRNA production)}\\
    &- \mathbf{\Gamma}_m \mathbf{\left\langle{m^{(x + 1)} p^y}\right\rangle} + 
    \mathbf{\Gamma}_m \mathbf{\left\langle{m (m - 1)^x p^y}\right\rangle}
    \text{  (mRNA degradation)}\\
    &- \mathbf{R}_p \mathbf{\left\langle{m^{(x + 1)} p^y}\right\rangle} + 
    \mathbf{R}_p \mathbf{\left\langle{m^{(x + 1)} (p + 1)^y}\right\rangle}
    \text{  (protein production)}\\
    &- \mathbf{\Gamma}_p \mathbf{\left\langle{m^x p^{(y + 1)}}\right\rangle} + 
    \mathbf{\Gamma}_p \mathbf{\left\langle{m^x p (p - 1)^y}\right\rangle}
    \text{  (protein degradation)}.
  \end{split}
\end{equation}

### Moment closure of the simple-repression distribution

A very interesting and useful feature of is that for a given value of $x$ and
$y$ the moment $\mathbf{\left\langle m^x p^y \right\rangle}$ is only a function
of lower moments. Specifically $\mathbf{\left\langle m^x p^y \right\rangle}$ is
a function of moments $\mathbf{\left\langle{m^{x'} p^{y'}}\right\rangle}$ that
satisfy two conditions: 
\begin{equation}
\begin{split}
    &1) y' \leq y,\\
  &2) x' + y' \leq x + y.
\end{split}
\end{equation}

To prove this we rewrite as 
\begin{equation}
\begin{split}
    \frac{d\mathbf{\left\langle m^x p^y \right\rangle}}{dt} &=
    \mathbf{K} \mathbf{\left\langle m^x p^y \right\rangle}\\
    &+ \mathbf{R}_m 
    \mathbf{\left\langle{p^y \left[ (m + 1)^x -m^x \right]}\right\rangle}\\
    &+ \mathbf{\Gamma}_m 
    \mathbf{\left\langle{m p^y \left[ (m - 1)^x - m^x \right]}\right\rangle}\\
    &+ \mathbf{R}_p 
    \mathbf{\left\langle{m^{(x + 1)}
    \left[ (p + 1)^y - p^y \right]}\right\rangle}\\
    &+ \mathbf{\Gamma}_p 
    \mathbf{\left\langle{m^x p \left[ (p - 1)^y - p^y \right]}\right\rangle},
\end{split}
\end{equation}
where the factorization is valid given the linearity of expected values. Now the
objective is to find the highest moment for each term once the relevant
binomial, such as $(m-1)^x$, is expanded. Take, for example, a simple case in
which we want to find the second moment of the mRNA distribution. We then set $x
= 2$ and $y = 0$. then becomes 
\begin{equation}
\begin{split}
    \frac{\mathbf{\left\langle{m^2 p^0}\right\rangle}}{dt} &=
    \mathbf{K} \mathbf{\left\langle{m^2 p^0}\right\rangle}\\
    &+ \mathbf{R}_m 
    \mathbf{\left\langle{p^0 \left[ (m + 1)^2 - m^2 \right]}\right\rangle}\\
    &+ \mathbf{\Gamma}_m
    \mathbf{\left\langle{m p^0 \left[ (m - 1)^2 - m^2 \right]}\right\rangle}\\
    &+ \mathbf{R}_p 
    \mathbf{\left\langle{m^{(2 + 1)} 
    \left[ (p + 1)^0 - p^0 \right]}\right\rangle}\\
    &+ \mathbf{\Gamma}_p 
    \mathbf{\left\langle{m^2 p \left[ (p - 1)^0 - p^0 \right]}\right\rangle}.
\end{split}
\end{equation}
Simplifying this equation gives 
\begin{equation}
\frac{d \mathbf{\left\langle{m^2}\right\rangle}}{dt} =
\mathbf{K} 
\mathbf{\left\langle{m^2}\right\rangle}
+ \mathbf{R}_m 
\mathbf{\left\langle{\left[ 2m + 1 \right]}\right\rangle}
+ \mathbf{\Gamma}_m 
\mathbf{\left\langle{\left[- 2m^2 + m \right]}\right\rangle}.
\end{equation}

XXX satisfies both of our conditions. Since we set $y$ to be zero, none of the
terms depend on any moment that involves the protein number, therefore $y' \leq
y$ is satisfied. Also the highest moment in also satisfies $x' + y' \leq x + y$
since the second moment of mRNA doesn't depend on any moment higher than
$\mathbf{\left\langle{m^2}\right\rangle}$. To demonstrate that this is true for
any $x$ and $y$ we now rewrite , making use of the binomial expansion
\begin{equation}
(z \pm 1)^n = \sum_{k=0}^n {n \choose k} (\pm 1)^{k} z^{n-k}.
\end{equation}
 Just as before let's look at each term individually. For the mRNA production
term we have 
\begin{equation}
\mathbf{R}_m 
\mathbf{\left\langle{p^y \left[ (m + 1)^x -m^x \right]}\right\rangle} =
\mathbf{R}_m 
\mathbf{\left\langle{p^y 
\left[ \sum_{k=0}^x {x \choose k} m^{x-k} - m^x \right]}\right\rangle}.
\end{equation}
When $k = 0$, the term inside the sum on the right-hand side cancels with the
other $m^x$, so we can simplify to
\begin{equation}
\mathbf{R}_m 
\mathbf{\left\langle{p^y \left[ (m + 1)^x -m^x \right]}\right\rangle} =
\mathbf{R}_m 
\mathbf{\left\langle{p^y 
\left[ \sum_{k=1}^x {x \choose k} m^{x-k} \right]}\right\rangle}.
\end{equation}
Once the sum is expanded we can see that the highest moment in this sum is given
by $\mathbf{\left\langle{m^{(x-1)} p^y}\right\rangle}$ which satisfies both of
the conditions on XXX.

For the mRNA degradation term we similarly have
\begin{equation}
\mathbf{\Gamma}_m 
\mathbf{\left\langle{m p^y \left[ (m - 1)^x - m^x \right]}\right\rangle} =
\mathbf{\Gamma}_m 
\mathbf{\left\langle{m p^y 
\left[ \sum_{k=0}^x {x \choose k}(-1)^k m^{x-k} - m^x \right]}\right\rangle}.
\end{equation}
Simplifying terms we obtain
\begin{equation}
\mathbf{\Gamma}_m 
\mathbf{\left\langle{m p^y \left[ \sum_{k=0}^x {x \choose k}(-1)^k m^{x-k} -
m^x \right]}\right\rangle} =
\mathbf{\Gamma}_m 
\mathbf{\left\langle{p^y 
\left[ \sum_{k=1}^x {x \choose k}(-1)^k m^{x+1-k} \right]}\right\rangle}.
\end{equation}
The largest moment in this case is $\mathbf{\left\langle{m^x p^y
}\right\rangle}$, which again satisfies the conditions on XXX.

The protein production term gives
\begin{equation}
\mathbf{R}_p 
\mathbf{\left\langle{m^{(x + 1)} \left[ (p + 1)^y - p^y \right]}\right\rangle} =
\mathbf{R}_p \mathbf{\left\langle{m^{(x + 1)}
\left[ \sum_{k=0}^y {y \choose k} (-1)^k p^{y-k} - p^y \right]}\right\rangle}.
\end{equation}
Upon simplification we obtain
\begin{equation}
\mathbf{R}_p 
\mathbf{\left\langle{m^{(x + 1)} 
\left[ \sum_{k=0}^y {y \choose k} (-1)^k p^{y-k} - p^y \right]}\right\rangle} =
\mathbf{R}_p 
\mathbf{\left\langle{m^{(x + 1)} 
\left[ \sum_{k=1}^y {y \choose k}(-1)^k p^{y-k} \right]}\right\rangle}.
\end{equation}
Here the largest moment is given by $\mathbf{\left\langle{m^{x+1}
p^{y-1}}\right\rangle}$, that again satisfies both of our conditions. For the
last term, for protein degradation we have
\begin{equation}
\mathbf{R}_p 
\mathbf{\left\langle{m^{(x + 1)} \left[ (p + 1)^y - p^y \right]}\right\rangle} =
\mathbf{R}_p 
\mathbf{\left\langle{m^{(x + 1)} \left[ \sum_{k=1}^y {y \choose k} (-1^k) p^{y - k}
  \right]}\right\rangle}.
\end{equation}
The largest moment involved in this term is therefore $\mathbf{\left\langle{m^x
p^{y-1}}\right\rangle}$. With this, we show that the four terms involved in our
general moment equation depend only on lower moments that satisfy XXX.

As a reminder, what we showed in this section is that the kinetic model
introduced in (A) has no moment-closure problem. In other words, moments of the
joint mRNA and protein distribution can be computed just from knowledge of lower
moments. This allows us to cleanly integrate the moments of the distribution
dynamics as cells progress through the cell cycle.

### Computing single promoter steady-state moments

(Note: The Python code used for the calculations presented in this section can
be found in the [following
link](https://www.rpgroup.caltech.edu//chann_cap/software/chemical_master_steady_state_moments_general.html)
as an annotated Jupyter notebook)

As discussed in XXX, one of the main factors contributing to cell-to-cell
variability in gene expression is the change in gene copy number during the cell
cycle as cells replicate their genome before cell division. Our minimal model
accounts for this variability by considering the time trajectory of the
distribution moments as given by . These predictions will be contrasted with the
predictions from a kinetic model that doesn't account for changes in gene copy
number during the cell cycle in XXX.

If we do not account for change in gene copy number during the cell cycle or for
the partition of proteins during division, the dynamics of the moments of the
distribution described in this section will reach a steady state. In order to
compute the steady-state moments of the kinetic model with a single gene across
the cell cycle, we use the moment closure property of our master equation. By
equating to zero for a given $\mathbf{x}$ and $\mathbf{y}$, we can solve the
resulting linear system and obtain a solution for $\mathbf{\left\langle m^x p^y
\right\rangle}$ at steady state as a function of moments
$\mathbf{\left\langle{m^{x'} p^{y'}}\right\rangle}$ that satisfy . Then, by
solving for the zero$^\text{th}$ moment $\mathbf{\left\langle{m^0
p^0}\right\rangle}$ subject to the constraint that the probability of the
promoter being in any state should add up to one, we can substitute back all of
the solutions in terms of moments $\mathbf{\left\langle{m^{x'}
p^{y'}}\right\rangle}$ with solutions in terms of the rates shown in . In other
words, through an iterative process, we can get at the value of any moment of
the distribution. We start by solving for the zero$^\text{th}$ moment. Since all
higher moments, depend on lower moments we can use the solution of the
zero$^\text{th}$ moment to compute the first mRNA moment. This solution is then
used for higher moments in a hierarchical iterative process.