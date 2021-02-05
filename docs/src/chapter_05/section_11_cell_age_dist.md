## Derivation of the cell age distribution 

E. O. Powell first derive in 1956 the distribution of cell age for a cell
population growing steadily in the exponential phase [@Powell1956]. This
distribution is of the form
\begin{equation}
P(a) = \ln(2) \cdot 2^{1 - a},
\end{equation}
where $a \in [0, 1]$ is the fraction of the cell cycle, 0 being the moment right
after the mother cell divides, and 1 being the end of the cell cycle just before
cell division. In this section we will reproduce and expand the details on each
of the steps of the derivation.

For an exponentially growing bacterial culture, the cells satisfy the growth law
\begin{equation}
{\frac{dn}{dt}} = \mu n,
\end{equation}
where $n$ is the number of cells and $\mu$ is the growth rate in units of
time$^{-1}$. We begin by defining $P(a)$ to be the probability density function
of a cell having age $a$. At time zero of a culture in exponential growth, i.e.
the time when we start considering the growth, not the initial condition of the
culture, there are $NP(a)da$ cells with age range between $[a, a + da]$. In
other words, for $N \gg 1$ and $da \ll a$
\begin{equation}
N P(a \leq x \leq a + da) \approx N P(a)da.
\end{equation}
We now define
\begin{equation}
F(\tau) = \int_\tau^\infty f(\xi) d\xi,
\end{equation}
as the fraction of cells whose division time is greater than $\tau$. This is
because in principle not all cells divide exactly after $\tau$ minutes, but
there is a distribution function $f(\tau)$ for the division time after birth.
Empirically it has been observed that a generalize Gamma distribution fits well
to experimental data on cell division time, but we will worry about this
specific point later on.

From the definition of $F(\tau)$ we can see that if a cell reaches an age $a$,
the probability of surviving to an age $a + t$ without dividing is given by
\begin{equation}
P(\text{age} = (a + t) \mid \text{age} = a) = F(a + t \mid a) =
\frac{F(a + t)}{F(a)}.
\end{equation}
This result comes simply from the definition of conditional probability. Since
$F(a)$ is the probability of surviving $a$ or more minutes without dividing, by
the definition of conditional probability we have that
\begin{equation}
F(a + t \mid a) = \frac{F(a, a + t)}{F(a)},
\end{equation}
where $F(a, a + t)$ is the joint probability of surviving $a$ minutes and $a +
t$ minutes. But the probability of surviving $a + t$ minutes or more implies
that the cell already survived $a$ minutes, therefore the information is
redundant and we have 
\begin{equation}
F(a, a + t) = F(a + t).
\end{equation}
This explains XXX. From this equation we can find that out of the $N P(a)da$
cells with age $a$ only a fraction
\begin{equation}
\left[ NP(a)da \right] F(a + t \mid a) = NP(a) \frac{F(a + t)}{F(a)} da
\end{equation}
will survive without dividing until time $a + t$. During that time interval $t$
the culture has passed from $N$ cells to $N e^{\mu t}$ cells given the
assumption that they are growing exponentially. The survivors $NP(a)F(a + t \mid
a)da$ then represent a fraction of the total number of cells
\begin{equation}
\frac{\text{\# survivors}}{\text{\# total cells}} =
\frac{\left[ NP(a)da \right] F(a + t \mid a)}{Ne^{\mu t}} =
  P(a)\frac{F(a + t)}{F(a)}da \frac{1}{e^{\mu t}},
\end{equation}
and their ages lie in the range $[a+t, a+t+da]$. Since we assume that the
culture is in steady state then it follows that the fraction of cells that
transitioned from age $a$ to age $a + t$ must be $P(a + t)da$. Therefore we have
a difference equation - the discrete analogous of a differential equation - of
the form
\begin{equation}
P(a + t) da = P(a) \frac{F(a + t)}{F(a)}e^{-\mu t} da.
\end{equation}
What this equation shows is a relationship that connects the probability of
having a life time of $a + t$ with a probability of having a shorter life time
$a$ and the growth of the population. If we take $t$ to be very small,
specifically if we assume $t \ll \mu^{-1}$ we can Taylor expand around $a$ the
following terms:
\begin{equation}
F(a + t) \approx F(a) + \frac{dF}{da} t,
\end{equation}
\begin{equation}
P(a + t) \approx P(a) + \frac{dP}{da} t,
\end{equation}
and
\begin{equation}
e^{-\mu t} \approx 1 - \mu t.
\end{equation}
Substituting these equations into gives
\begin{equation}
P(a) + \frac{dP}{da} t = P(a) \left( \frac{F(a) + \frac{dF}{da}t}{
  F(a)} \right) (1 - \mu t).
\end{equation}
This can be rewritten as
\begin{equation}
\frac{1}{P(a)} \frac{dP}{da} =
\frac{1}{F(a)} \frac{dF}{da} - \mu - \frac{\mu t}{F(a)} \frac{dF}{da}.
\end{equation}
Since we assumed $t \ll \mu^{-1}$ we then approximate the last term to be close
to zero. We can then simplify this result into
\begin{equation}
\frac{1}{P(a)} \frac{dP}{da} = \frac{1}{F(a)} \frac{dF}{da} - \mu.
\end{equation}
Integrating both sides of the equation with respect to $a$ gives
\begin{equation}
\ln P(a) = \ln F(a) - \mu a + C,
\end{equation}
where $C$ is the integration constant. Exponentiating both sides gives 
\begin{equation}
P(a) = C' F(a)e^{-\mu a}.
\end{equation}
Where $C' \equiv e^C$. To obtain the value of the unknown constant we recall
that $F(0) = 1$ since the probability of having a life equal or longer than zero
must add up to one, therefore we have that $P(0) = C'$. This gives then 
\begin{equation}
P(a) = P(0) e^{-\mu a} F(a).
\end{equation}
Substituting the definition of $F(a)$ gives
\begin{equation}
P(a) = P(0) e^{-\mu a} \int_a^\infty f(\xi) d\xi.
\end{equation}
The last step of the derivation involves writing $P(0)$ and the growth rate
$\mu$ in terms of the cell cycle length distribution $f(\tau)$.

The growth rate of the population cell number (not the growth of cell mass) is
defined as the number of cell doublings per unit time divided by the number of
cells. This is more clear to see if we write as a finite difference
\begin{equation}
\frac{N(t + \Delta t) - N(t)}{\Delta t} = \mu N(t).
\end{equation}
If the time $\Delta t$ is the time interval it takes to go from $N$ to $2N$
cells we have 
\begin{equation}
\frac{2N - N}{\Delta t} = \mu N.
\end{equation}
Solving for $\mu$ gives
\begin{equation}
\mu = \overbrace{\frac{2N - N}{\Delta t}}
^{\text{\# doubling events per unit time}}
\overbrace{\frac{1}{N}}^{\frac{1}{\text{population size}}}.
\end{equation}
We defined $F(a)$ to be the probability of a cell reaching an age $a$ or
greater. For a cell to reach an age $a + da$ we can then write
\begin{equation}
F(a + da) = \int_{a + da}^{\infty} f(\xi) d\xi
= \int_a^{\infty} f(\xi) d\xi - \int_a^{a + da} f(\xi) d\xi.
\end{equation}
We can approximate the second term on the right hand side to be
\begin{equation}
\int_a^{a + da} f(\xi) d\xi \approx f(a) da,
\end{equation}
for $da \ll a$, obtaining 
\begin{equation}
F(a + da) \approx F(a) - f(a)da.
\end{equation}
What this means is that from the original fraction of cells $F(a)$ with age $a$
or greater a fraction $f(a)da / F(a)$ will not reach age $(a + da)$ because they
will divide. So out of the $NP(a)$ cells that reached exactly age $a$, the
number of doubling events on a time interval $da$ is given by
\begin{equation}
{\text{\# doublings of cells of age } a {\text{ on interval } da}} =
  \overbrace{NP(a)}^{\text{\# cells of age }a}
  \overbrace{\frac{f(a) da}{F(a)}}^{\text{fraction of doublings per unit time}}.
\end{equation}
The growth rate then is just the sum (integral) of each age contribution
to the total number of doublings. This is
\begin{equation}
\mu = \frac{1}{N} \int_0^\infty NP(a) \frac{f(a)da}{F(a)}.
\end{equation}
Substituting gives
\begin{equation}
\mu = \int_0^\infty [P(0) e^{-\mu a} F(a)] \frac{f(a)da}{F(a)}
  = \int_0^\infty P(0) e^{-\mu a} f(a)da.
\end{equation}
We now have the growth rate $\mu$ written in terms of the cell cycle length
probability distribution $f(a)$ and the probability $P(0)$. Since $P(a)$ is a
probability distribution it must be normalized, i.e. 
\begin{equation}
\int_0^\infty P(a) da = 1.
\end{equation}
Substituting into this normalization constraint gives
\begin{equation}
\int_0^\infty P(0) e^{-\mu a} F(a) da = 1.
\end{equation}
From here we can
integrate the left hand side by parts. We note that given the definition
of $F(a)$, the derivative with respect to $a$ is $-f(a)$ rather than
$f(a)$. This is because if we write the derivative of $F(a)$ we have
\begin{equation}
\frac{dF(a)}{da} \equiv \lim_{da \rightarrow 0}
  \frac{F(a + da) - F(a)}{da}.
\end{equation}
Substituting the definition of $F(a)$ gives 
\begin{equation}
\frac{dF(a)}{da} = \lim_{da \rightarrow 0} \frac{1}{da}
\left[\int_{a + da}^\infty f(\xi) d\xi - \int_a^\infty f(\xi) d\xi \right].
\end{equation}
This difference in the integrals can be simplified to
\begin{equation}
\lim_{da \rightarrow 0} \frac{1}{da} \left[ \int_{a + da}^\infty f(\xi) d\xi -
  \int_a^\infty f(\xi) d\xi \right]\approx \frac{-f(a)da}{da} = -f(a).
\end{equation}
Taking this into account we now perform the integration by parts obtaining 
\begin{equation}
P(0) \left[ \frac{e^{-\mu t}}{-\mu} F(a) \right]^\infty_0
 - P(0) \int_0^\infty \frac{e^{-\mu a}}{-\mu} (-f(a)) da = 1.
\end{equation}
On the first term on the left hand side we have that as $a \rightarrow \infty$,
both terms $e^{-\mu a}$ and $F(a)$ go to zero. We also have that $e^{\mu 0} = 1$
and $F(0) = 1$. This results in
\begin{equation}
\frac{P(0)}{\mu} - P(0) \int_0^\infty \frac{e^{-\mu a}}{\mu} f(a) da = 1.
\end{equation}
The second term on the left hand side is equal to since
\begin{equation}
\mu = \int_0^\infty P(0) e^{-\mu a} f(a)da \Rightarrow
  1 = \int_0^\infty P(0) \frac{e^{-\mu a}}{\mu} f(a)da.
\end{equation}
This implies that on we have 
\begin{equation}
\frac{P(0)}{\mu} - 1 = 1 \Rightarrow P(0) = 2 \mu.
\end{equation}
With this result in hand we can rewrite as
\begin{equation}
P(a) = 2\mu e^{-\mu a} \int_a^\infty f(\xi) d\xi.
\end{equation}
Also we can rewrite the result for the growth rate $\mu$ on as
\begin{equation}
\mu = 2 \mu \int_0^\infty e^{-\mu a} f(a) da \Rightarrow
  2 \int_0^\infty e^{-\mu a} f(a) da = 1.
\end{equation}

As mentioned before the distribution $f(a)$ has been empirically fit to a
generalize Gamma distribution. But if we assume that our distribution has almost
negligible dispersion around the mean average doubling time $a = \tau_d$, we can
approximate $f(a)$ as
\begin{equation}
f(a) = \delta(a - \tau_d),
\end{equation}
a Dirac delta function. Applying this to results in 
\begin{equation}
2 \int_0^\infty e^{-\mu a} \delta(a - \tau_a) da = 1
  \Rightarrow 2 e^{-\mu \tau_d} = 1.
\end{equation}
Solving for $\mu$ gives
\begin{equation}
\mu = \frac{\ln 2}{\tau_d}.
\end{equation}
This delta function approximation for $f(a)$ has as a consequence that 
\begin{equation}
F(a) =
  \begin{cases}
    1 \text{ for } a \in [0, \tau_d],\\
    0 \text{ for } a > \tau_d.
  \end{cases}
\end{equation}
Fianlly we can rewrite as
\begin{equation}
P(a) = 2 \left( \frac{\ln 2}{\tau_d} \right)
e^{- \frac{\ln 2}{\tau_d} a} \int_a^\infty \delta(\xi - \tau_d) d\xi
\Rightarrow = 2 \ln 2 \cdot 2^\frac{-a}{\tau_d}.
\end{equation}
Simplifying this we obtain 
\begin{equation}
P(a) =
  \begin{cases}
    \ln 2 \cdot 2^{1 - \frac{a}{\tau_d}} \text{ for } a \in [0, \tau_d],\\
    0 \text{ otherwise}.
  \end{cases}
\end{equation}
This is the equation we aimed to derive. The distribution of cell ages over the
cell cycle.
