### Computing the moments of the mRNA and protein distributions 

Finding analytical solutions to chemical master equations is often fraught with
difficulty. An alternative approach is to to approximate the distribution. One
such scheme of approximation, the maximum entropy principle, makes use of the
moments of the distribution to approximate the full distribution. In this
section we will demonstrate an iterative algorithm to compute the mRNA and
protein distribution moments.

The kinetic model for the simple repression motif depicted in
[@Fig:ch3_fig02](A) consists of an infinite system of ODEs for each possible
pair of mRNA and protein copy number, $(m, p)$. To compute any moment of the
distribution, we define a vector
$$
\langle \mathbf{m^x p^y} \rangle \equiv 
(\langle m^x p^y\rangle_A, 
\langle m^x p^y \rangle_I, 
\langle m^x p^y\rangle_R)^T,
\label{eq:ch2_eq07}
$$
where $\langle m^x p^y \rangle_S$ is the expected value of $m^x p^y$ in state $S
\in \{A, I, R\}$ for $x, y \in \mathbb{N}$. In other words, just as we defined
the vector $\mathbf{P}(m, p)$, here we define a vector to collect the expected
value of each of the promoter states. By definition, any of these moments
$\langle m^x p^y \rangle_S$ can be computed as
$$
\langle m^x p^y \rangle_S \equiv 
\sum_{m=0}^\infty \sum_{p=0}^\infty m^x p^y P_S(m, p).
\label{eq:ch3_eq08}
$$
Summing over all possible values for $m$ and $p$ in Eq. $\ref{eq:ch3_eq06}$
results in an ODE for any moment of the distribution of the form (See Chapter 5
for full derivation) 
$$
\begin{split}
    \frac{d \mathbf{\langle m^x p^y \rangle }}{dt} &=
    \mathbf{K} \mathbf{\langle m^x p^y \rangle}\\
    &+ \mathbf{R}_m \mathbf{\langle p^y \left[ (m + 1)^x -m^x \right]\rangle }
     + \mathbf{\Gamma}_m \mathbf{\langle m p^y \left[ (m - 1)^x - m^x \right]\rangle}\\
    &+ \mathbf{R}_p \mathbf{\langle m^{(x + 1)} \left[ (p + 1)^y - p^y \right]\rangle}
     + \mathbf{\Gamma}_p \mathbf{\langle m^x p \left[ (p - 1)^y - p^y \right]\rangle }.
\end{split}
\label{eq:ch3_eq09}
$$

Given that all transitions in our stochastic model are first order reactions,
has no moment-closure problem [@Voliotis2014a]. This means that the dynamical
equation for a given moment only depends on lower moments (See Chapter 5 for
full proof). This feature of our model implies, for example, that the second
moment of the protein distribution $\langle p^2 \rangle$ depends only on the
first two moments of the mRNA distribution $\langle m \rangle$ and $\langle m^2
\rangle$, the first protein moment $\langle p \rangle$, and the
cross-correlation term $\langle mp \rangle$. We can therefore define
$\boldsymbol{\mu}^{\mathbf{(x, y)}}$ to be a vector containing all moments up to
$\mathbf{\langle m^x p^y\rangle}$ for all promoter states, 
$$
\boldsymbol{\mu}^{\mathbf{(x, y)}} = \left[ \mathbf{\langle m^0 p^0 \rangle},
\mathbf{\langle m^1 p^0 \rangle},
\ldots, \mathbf{\langle m^x p^y \rangle} \right]^T.
\label{eq:ch3_eq10}
$$
Explicitly for the three-state promoter model depicted in [@Fig:ch3_fig02](A)
this vector takes the form
$$
\boldsymbol{\mu}^{\mathbf{(x, y)}} = 
\left[ 
    \langle m^0 p^0 \rangle_A,
    \langle m^0 p^0 \rangle_I,
    \langle m^0 p^0 \rangle_R,
    \ldots,
    \langle m^x p^y \rangle_A,
    \langle m^x p^y \rangle_I,
    \langle m^x p^y \rangle_R 
\right]^T.
\label{eq:ch3_eq11}
$$

Given this definition we can compute the general moment dynamics as 
$$
\frac{d \boldsymbol{\mu}^{\mathbf{(x, y)}}}{dt} = \mathbf{A}
\boldsymbol{\mu}^{\mathbf{(x, y)}}, 
\label{eq:ch3_eq12}
$$
where $\mathbf{A}$ is a square matrix that contains all the numerical
coefficients that relate each of the moments. We can then use Eq. $\ref{eq:ch3_eq09}$ to
build matrix $\mathbf{A}$ by iteratively substituting values for the exponents
$x$ and $y$ up to a specified value. In the next section, we will use
Eq. $\ref{eq:ch3_eq12}$ to numerically integrate the dynamical equations for our moments
of interest as cells progress through the cell cycle. We will then use the value
of the moments of the distribution to approximate the full gene expression
distribution. This method is computationally more efficient than trying to
numerically integrate the infinite set of equations describing the full
probability distribution $\mathbf{P}(m, p)$, or using a stochastic algorithm to
sample from the distribution.
