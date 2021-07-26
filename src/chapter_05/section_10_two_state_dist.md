## Derivation of the steady-state mRNA distribution

In this section, we will derive the two-state promoter mRNA distribution we
quote in [Sec. 5.2](#sec:ch5_sec03). For this method, we will make use of the
so-called generating functions. Generating functions are mathematical objects on
which we can encode a series of infinite numbers as coefficients of a power
series. The power of generating functions comes from the fact that we can
convert an infinite-dimensional system of coupled ordinary differential
equations--in our case, the system of differential equations defining all
probabilities $P(m, t)$ for $m \in \mathbb{Z}$--into a single partial
differential equation that we can then solve to extract back the probability
distributions.

To motivate the use of generating functions, we will begin with the simplest
case: the one-state Poisson promoter.

### One-state Poisson promoter

We begin by defining the reaction scheme that defines the one-state promoter.
[@Fig:ch5_fig35] shows the schematic representation of the Poisson promoter as a
simple cartoon (part (A)) and as the Markov chain that defines the state space
of the system (part (B)).

![**One-state Poisson promoter.** (A) Schematic of the kinetics of the one
state-promoter. mRNA is produced and degrade stochastically with a rate $r_m$
and $\gamma_m$, respectively. (B) Representation of the Markov-chain for the
state space that the promoter can be. The distribution $P(m, t)$ represents the
probability of having a certain discrete number of mRNA $m$ at time $t$. The
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
\label{eq:one_state_master}
$$

When solving for the distribution, our objective is to obtain the equation that
defines $P(m, t)$ for all possible values of $m \in \mathbb{Z}$. The power of
the generating functions is that these probability distribution values are used
as a power series's coefficients. To make this clear, let us define the
generating function $G(z, t)$ as
$$
G(z, t) \equiv \sum_{m=0}^\infty z^m P(m, t),
\label{eq:gen_fn_def}
$$
where $z$ is a "dummy" variable that we don't care about. The reason this is
useful is that if we find the closed-form solution for this generating function,
and we are able to split the factor $z^m$ from its coefficient $P(m, t)$, then
we will have find the solution for the distribution. Furthermore, the generating
function allows us to compute the moments of the distribution. For example, for
the zeroth moment $\langle m^0 \rangle$ we know that
$$
\langle m^0 \rangle = \sum_{m=0}^\infty m^0 P(m, t) = 1,
$$
i.e., this is the normalization constraint of the distribution. From the
definition of the generating function we can then see that
$$
G(1, t) = \sum_{m=0}^\infty 1^m P(m, t) = 1.
\label{eq:generating_norm}
$$
Furthermore, the first moment of the distribution is defined as
$$
\langle m \rangle = \sum_{m=0}^\infty m P(m, t).
$$
From the definition of the generating function, we can construct this quantity
by computing
$$
\left. \frac{\partial G(z, t)}{\partial z} \right\vert_{z=1} =
\frac{\partial}{\partial z} \left[ 
\sum_{m=0}^\infty z^m P(m, t)
\right]_{z=1} =
\sum_{m=0}^\infty m P(m, t).
$$
Therefore we have that
$$
\langle m \rangle = 
\left. \frac{\partial G(z, t)}{\partial z} \right\vert_{z=1}.
\label{eq:first_mom_gen}
$$
Similar constructions can be built for higher moments of the distribution.

Let us then apply the definition of the generating function to Eq.
$\ref{eq:one_state_master}$. For this, we multiply both sides by $z^m$ and sum
over all values of $m$, obtaining
$$
\begin{split}
\sum_{m=0}^\infty z^m
\frac{d P(m, t)}{dt} &=
\sum_{m=0}^\infty z^m \left[
r_m P(m-1, t)
- r_m P(m, t) \right. \\
&+ \left. \gamma_m (m+1) P(m+1, t)
- \gamma_m m P(m, t).
\right]
\end{split}
$$
Distributing the sum, we find
$$
\begin{split}
\frac{d}{dt} \sum_{m=0}^\infty z^m P(m, t) &=
\sum_{m=0}^\infty z^m r_m P(m-1, t)
- \sum_{m=0}^\infty z^m r_m P(m, t) \\
&+ \sum_{m=0}^\infty z^m \gamma_m (m+1) P(m+1, t)
- \sum_{m=0}^\infty z^m \gamma_m m P(m, t).
\end{split}
\label{eq:one_state_master_sum}
$$
We see that the terms involving $z^m P(m, t)$ can be directly substituted with
Eq. $\ref{eq:gen_fn_def}$. For the other terms, we have to be slightly more
clever. The first trick will allow us to rewrite the term involving $z^m m P(m,
t)$ as
$$
\begin{aligned}
\sum_{m} z^{m} \cdot m \cdot P(m, t) &=
\sum_{m} z \frac{\partial z^{m}}{\partial z} P(m, t), \\
&=\sum_{m} z \frac{\partial}{\partial z}\left(z^{m} P(m, t)\right), \\
&=z \frac{\partial}{\partial z}\left(\sum_{m} z^{m} P(m, t)\right), \\
&=z \frac{\partial G(z, t)}{\partial z}.
\end{aligned}
$$
Next, let us deal with the term involving $(m+1)$. We first define $k = m + 1$.
With this, we can write
$$
\begin{aligned}
\sum_{m=0}^{\infty} z^{m} \cdot(m+1) \cdot P(m+1, t) &=
\sum_{k=1}^{\infty} z^{k-l} \cdot k \cdot P(k, t), \\
&=z^{-1} \sum_{k=1}^{\infty} z^{k} \cdot k \cdot P(k, t), \\
&=z^{-1} \sum_{k=0}^{\infty} z^{k} \cdot k \cdot P(k, t), \\
&=z^{-1}\left(z \frac{ \partial G(z)}{\partial z}\right), \\
&=\frac{\partial G(z)}{\partial z},
\end{aligned}
$$
where for the third step, we reindexed the sum to include $k=0$ since it does
not contribute to the total sum. Finally, for the term involving $P(m-1, t)$.
For this we define $k = m-1$. This allows us to rewrite the term as
$$
\begin{aligned}
\sum_{m=0}^{\infty} z^{m} P(m-1, t) &=\sum_{k=-1}^{\infty} z^{k+1} P(k, t), \\
&=\sum_{k=0}^{\infty} z^{k+1} P(k, t), \\
&=z \sum_{k=0}^{\infty} z^{k} P(k, t), \\
&=z G(z, t)
\end{aligned}
$$
For the second step we reindexed the sum from $-1$ to $0$ since $P(-1, t) = 0$.

All of these clever reindexing allows us to rewrite Eq.
$\ref{eq:one_state_master_sum}$ as
$$
\frac{\partial G(z, t)}{\partial t} =
r z G(z, t)-r G(z, t)
+ \gamma \frac{\partial G(z, t)}{\partial z}
- \gamma z \frac{\partial G(z, t)}{\partial z}.
$$
Factorizing terms we have
$$
\frac{\partial G(z, t)}{\partial t} 
=-r G(z, t)(1-z)
+\gamma \frac{\partial G(z, t)}{\partial z}(1-z).
$$
Let us appreciate how beautiful this is: we took an infinite-dimensional system
of ordinary differential equations--the master equation--and turn it into a
single partial differential equation (PDE). All we have to do now is solve this
PDE, and then transform the solution into a power series to extract the
distribution.

Let us focus on the steady-state case. For this, we set the time derivative to
zero. Doing this cancels the $(1-z)$ term, leaving a straightforward ordinary
differential equation for $G(z)$
$$
\frac{dG(z)}{dz} = \frac{r}{\gamma} G(z).
$$
Solving this equation by separation of variables results in
$$
G(z) = C e^{\frac{r}{\gamma}z}.
$$
To obtain the integration constant, we use the normalization condition of the
probability distribution (Eq. $\ref{eq:generating_norm}$), obtaining
$$
1 = C e^{\frac{r}{\gamma}} \Rightarrow
C = e^{-\frac{r}{\gamma}}.
$$
This means that the generating function takes the form
$$
G(z) = e^{-\frac{r}{\gamma}} e^{\frac{r}{\gamma}z}.
$$
All we have left is trying to rewrite the generating function as a power series
on $z$. If we succeed in doing so, we will have recovered the probability
distribution $P(m, t)$. For this, we simply use the Taylor expansion of $e^x$,
obtaining
$$
G(z) = e^{-\frac{r}{\gamma}} 
\sum_{m=0}^\infty \frac{\left( \frac{r}{\gamma}z \right)^m}{m!}.
$$
From this form, it becomes clear how to split the $z^m$ term from the
coefficient that, by the definition of the generating function, is the
probability distribution we are looking for. The separation takes the form
$$
G(z)=\sum_{m=0}^{\infty} z^{m}
\left[\frac{e^{-r / \gamma}\left(\frac{r}{\gamma}\right)^{m}}{m !}\right],
$$
where we can see that we recover the expected Poisson distribution for this 
one-state promoter
$$
P(m) = e^{-r / \gamma}\frac{\left(\frac{r}{\gamma}\right)^{m}}{m !}.
$$

### Two-state promoter

Having shown the generating function's power, let us now turn our attention to
the relevant equation we are after: the two-state mRNA distribution. This model
assumes that the promoter can exist in two discrete states (See
[@Fig:ch5_fig36](A)): a transcriptionally active state $A$ from which
transcription can take place at a constant rate $r_m$, and an inactive state $I$
where no transcription takes place. The mRNA is stochastically degraded with a
rate $\gamma_m$ regardless of the state of the promoter. [@Fig:ch5_fig36](B)
shows the Markov chain that connects all of the possible states of the promoter.
For this particular case, there are not only "horizontal" transitions where the
mRNA copy number changes, but "vertical" transitions where only the promoter's
state changes. Because of this, we need to define two coupled master equations
that take the form
$$
\begin{aligned}
\frac{d P_{A}(m, t)}{d t} &=-k^{(p)}_{\text{off}} P_{A}(m, t) +
k^{(p)}_{\text{on}} P_{I}(m, t) \\
&+\gamma_m (m+1) P_{A}(m+1, t) - \gamma_m m P_{n}(m, t) \\
&+r_m P_{A}(m-1, t)-r_m P_{A}(m, t) \\
\end{aligned}
$$
for the active state, and
$$
\begin{aligned}
\frac{d P_{I}(m, t)}{dt} &=k^{(p)}_{\text{off}} P_{A}(m, t)-
k^{(p)}_{\text{on}} P_{I}(m, t) \\
&+\gamma_m (m+1) P_{I}(m+1, t)-\gamma_m m P_{I}(m, t),
\end{aligned}
$$
for the inactive state.

![**Two-state Poisson promoter.** (A) Schematic of the kinetics of the two-state
promoter. The promoter is imagined to exist in two-state--a transcriptionally
active state $A$ and an inactive state $I$. The transition between these states
is governed by the rates $k^{(p)}_{\text{on}}$ and $k^{(p)}_{\text{off}}$ mRNA
is produced and degrade stochastically with a rate $r_m$ and $\gamma_m$,
respectively. (B) Representation of the Markov-chain for the state space that
the promoter can be in. The distribution $P(m, t)$ represents the probability of
having a certain discrete number of mRNA $m$ at time $t$. The transition between
states depends on the previously mentioned rates.](ch5_fig36){#fig:ch5_fig36
short-caption="One-state Poisson promoter"}

#### Obtaining the partial differential equation for the generating function

The first thing we must do is to transform this infinite-dimensional system of
ordinary differential equations in $m$ to a single partial differential equation
using the generating function. For this particular case, there are two
generating functions of the form
$$
G_x(z, t) = \sum_{m=0}^\infty z^m P_x(m, t),
$$
where $x \in \{A, I \}$. The probability of having $m$ mRNA at time $t$
regardless of the promoter state is given by
$$
P(m, t) = P_A(m, t) + P_I(m, t).
$$
Therefore, the corresponding generating function for the whole system is given
by
$$
G(z, t) = G_A(z, t) + G_I(z, t).
$$

As with the one-state promoter case, let us transform our master equations by
multiplying both sides by $z^m$ and sum over all $m$. For the active state $A$
we have
$$
\begin{aligned}
\sum_{m} z^{m} \frac{d P_{A}(m, t)}{d t} =
\sum_{m} z^{m} &\left[-k^{(p)}_{\text{off}} P_{A}(m, t) 
+ k^{(p)}_{\text{on}} P_{I}(m, t)\right.\\
&+ \gamma_m (m+1) P_{A}(m+1, t)-\gamma_m m P_{m}(m, t) \\
&\left. + r_m P_{A}(m-1, t) - r_m P_{A}(m, t)\right].
\end{aligned}
$$
After distributing the sum, we can use the tricks from the previous, allowing
us to write this as a partial differential equation of the form
$$
\begin{aligned}
\frac{\partial G_{A}(z, t)}{\partial t} &=
-k^{(p)}_{\text{off}} G_{A}(z, t)+k^{(p)}_{\text{on}} G_{I}(z, t) \\
&-\gamma_m(z-1) \frac{\partial G_A(z, t)}{\partial z} 
+r_m(z-1) G_{A}(z, t).
\end{aligned}
\label{eq:gn_fn_act}
$$
An equivalent process can be done for the inactive state I, obtaining
$$
\begin{aligned}
\frac{\partial G_{I}(z, t)}{\partial t} &=
k^{(p)}_{\text{off}} G_{A}(z, t) - k^{(p)}_{\text{on}} G_{I}(z, t) \\
&-\gamma_m(z-1) \frac{\partial G A(z, t)}{\partial z} 
+r_m(z-1) G_{I}(z, t).
\end{aligned}
\label{eq:gn_fn_inact}
$$
We turned the infinite-dimensional system of ordinary differential equations
into a system of two coupled partial differential equations. Let us transform
the equations further. Since we have a common term $(z - 1)$, it will be
convenient to define $v \equiv (z -1)$. From the chain rule, it follows that
$$
d v=d(z-1)=d z \Rightarrow \frac{\partial G}{\partial v} = 
\frac{\partial G}{\partial z} \frac{d z}{d v}.
$$
Making this substitution in Eqs. $\ref{eq:gn_fn_act}$ and $\ref{eq:gn_fn_inact}$
results in
$$
\begin{aligned}
\frac{\partial G_{A}(v, t)}{\partial t} &=-k^{(p)}_{\text{off}} G_{A}(v, t)
+ k^{(p)}_{\text{on}} G_{I}(v, t) \\
&-\gamma_m v \frac{\partial G_{A}(v, t)}{\partial v} + r_m v G_{A}(v, t) \\
\end{aligned}
$$
for the actives state, and
$$
\begin{aligned}
\frac{\partial G I(v, t)}{\partial t}=& k^{(p)}_{\text{off}} G_{A}(v, t)
- k^{(p)}_{\text{on}} G_{I}(v, t) \\
&-r_m v \frac{\partial G_{I}(v, t)}{\partial v},
\end{aligned}
$$
for the active state. 

Since we care about the steady-state distribution, it is at this point that we
set the time derivative of both equations to zero. Doing this results in
$$
\gamma_{m} v \frac{d G_{A}(v)}{d v}= 
-k^{(p)}_{\text{off}} G_{A}(v)
+ k^{(p)}_{\text{on}} G_{I}(v).
+ r_m v G_{A}(v),
\label{eq:steady_act}
$$
and
$$
\gamma_{m} v \frac{d G_{I}(v)}{d v}= 
k^{(p)}_{\text{off}} G_{A}(v)
- k^{(p)}_{\text{on}} G_{I}(v).
\label{eq:steady_inact}
$$
Adding Eqs. $\ref{eq:steady_act}$ and $\ref{eq:steady_inact}$ gives a simple
result
$$
\gamma_m \frac{d G(v)}{dv} = r_m G_A(v).
\label{eq:gen_fn_rel}
$$

Our objective is not to write Eqs. $\ref{eq:steady_act}$ and
$\ref{eq:steady_inact}$ as a function of only one of the generating functions,
i.e., we want two independent differential equations. These equations are both
function of $G_A(v)$ and $G_I(v)$, but Eq. $\ref{eq:gen_fn_rel}$ tells us how to
relate both generating functions via the first derivative. This suggests that
taking another derivative of Eqs. $\ref{eq:steady_act}$ and
$\ref{eq:steady_inact}$ with respect to $z$ could be useful. Let us go ahead and
compute these derivatives. For the active state, we find
$$
\small
\gamma_m \frac{d G_{A}(v)}{d v} 
+ \gamma_{m} v \frac{d^{2} G_{A}(v)}{d v^{2}} =
-k^{(p)}_{\text{off}} \frac{d G_{A}(v)}{d v}
+ k^{(p)}_{\text{on}} \frac{d \sigma I(v)}{d v}
+ r_m G_{A}(v, t) + r_m v \frac{d G_{A}(v)}{d v}
$$
Upon simplification, we can write this Eq. as
$$
\gamma_m v \frac{d^2 G_A}{d v^{2}}
+ \left(\gamma_m +k^{(p)}_{\text{off}}
-r_m v\right) \frac{d G_A}{d v}
-k^{(p)}_{\text{on}} \frac{d G_I}{d v}
-r_m G_{A}(v)=0.
\label{eq:gen_fn_2nd_act}
$$
From Eq. $\ref{eq:gen_fn_rel}$ we have that
$$
\frac{G_I}{dv} = \frac{r_m}{\gamma_m} G_A(v) 
- \frac{d G_A}{dv}.
$$
Substituting this into $\ref{eq:gen_fn_2nd_act}$ results in
$$
\gamma_m v \frac{d^{2} G_A}{d v^{2}}
+ \left(\gamma_m + k^{(p)}_{\text{off}} 
+ k^{(p)}_{\text{on}} 
- r_m v \right) \frac{d G_A}{d v}
- r_m \left(1+\frac{k^{(p)}_{\text{on}}}{\gamma_m}\right) G_A(v) = 0.
\label{eq:gen_2nd_act_final}
$$
For the inactive state, upon taking a derivative with respect to $v$, we find
$$
\gamma_{m} v \frac{d^{2} G_{I}}{d v^{2}}
+ \left(\gamma_m + k^{(p)}_{\text{on}}\right) \frac{d G_I}{d v}
- k^{(p)}_{\text{off}} \frac{d G_A}{d v} = 0.
\label{eq:gen_fn_2nd_inact}
$$
Again from $\ref{eq:gen_fn_rel}$ we have that
$$
\frac{d G_A}{dv} = \frac{r_m}{\gamma_m} G_A - \frac{d G_I}{dv}.
$$
Substituting this result into Eq. $\ref{eq:gen_fn_2nd_inact}$ gives
$$
\gamma_{m} v \frac{d^{2} G_{I}}{d v^{2}}
+ \left(\gamma_m+k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}\right) 
\frac{d G_I}{d v} 
- \frac{k^{(p)}_{\text{off}} r_m}{\gamma_m} G_{A}(v)=0.
\label{eq:gen_fn_2nd_inact_2}
$$
So far, we have not removed the dependence on $G_A(v)$. But we notice that from
Eq. $\ref{eq:steady_inact}$ we have that
$$
G_A(v) = \frac{\gamma_m v}{k^{(p)}_{\text{off}}} \frac{d G_I}{d v}
+ \frac{k^{(p)}_{\text{on}}}{k^{(p)}_{\text{off}}} G_I.
$$
Using this identity allows us to write Eq. $\ref{eq:gen_fn_2nd_inact_2}$ as
$$
\gamma_{m} v \frac{d^{2} G_{I}}{d v^{2}}
+ \left(\gamma_m + k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}} -r_m v\right) 
\frac{d G_I}{d v} - \frac{k^{(p)}_{\text{on}} r_m}{\gamma_m} G_I = 0.
\label{eq:gen_2nd_inact_final}
$$

To obtain a single partial differential equation we add Eqs.
$\ref{eq:gen_2nd_act_final}$ and $\ref{eq:gen_2nd_inact_final}$, obtaining
$$
\gamma_m v \frac{d^{2} G}{d v^{2}}
+ \left(\gamma_m + k^{(p)}_{\text{off}} + k^{(p)}_{\text{on}} -r_m v\right) 
\frac{d G}{d v} 
- \frac{r_m k^{(p)}_{\text{on}}}{\gamma_m} G(v)
- r_m G_A(v) = 0,
$$
where we substituted $G_A(v) + G_I(v) = G(v)$. To remove the last $G_A(v)$ we
utilize again Eq. $\ref{eq:gen_fn_rel}$, obtaining
$$
\gamma_m v \frac{d^2 G}{dv^2} 
+ \left( k^{(p)}_{\text{off}} + k^{(p)}_{\text{on}} - r_m v \right)
\frac{dG}{dv} 
- \frac{r_m k^{(p)}_{\text{on}}}{\gamma_m} G(v) = 0.
\label{eq:gen_2nd}
$$

#### Solving the partial differential equation

Eq. $\ref{eq:gen_2nd}$ looks almost like the so-called Kummer's equation also
known as the confluent hypergeometric differential equation--a second order
differential equation of the form
$$
z \frac{d^2w}{dz^2} + (b - z) \frac{dw}{dz} - aw = 0.
\label{eq:kummer}
$$
The solution fo the Kummer equation can be expressed as the sum of two 
functions:
1. The confluent hypergeometric function of the first kind,
2. The Tricomi function.
This is written as
$$
w(z) = A {}_1F_1(a, b, z) + B z^{1-b} {}_1 F_1(a+1-b, 2-b, z),
\label{eq:kummer_sol}
$$
where $A$ and $B$ are constants, and ${}_1F_1$ is the confluent hypergeometric
function of the first kind defined as
$$
{}_1F_1(a, b, z) = \sum_{m=0}^{\infty} 
\frac{a^{(m)}z^n}{b^{(m)} m!},
$$
where $a^{(n})$ and $b^{(n)}$ are the rising factorials, i.e.,
$$
a^{(0)} = 1,
$$
and
$$
a^{(n)} = a (a + 1) (a + 2) \cdots (a + n - 1).
$$

To write Eq. $\ref{eq:gen_2nd}$ in the form of Eq. $\ref{eq:kummer}$ we can
define $s \equiv r_m v / \gamma_m$. The chain rule tells us that
$$
ds = \frac{r_m}{\gamma_m} dv \Rightarrow 
\frac{dG}{ds} = \frac{dG}{dv}\frac{dv}{ds} = 
\frac{\gamma_m}{r_m} \frac{dG}{dv}.
$$
From the chain rule, we also conclude that
$$
\frac{d^2G}{ds^2} = 
\frac{d}{dv} \left( \frac{dG}{dv} \frac{dv}{ds} \right) \frac{dv}{ds} =
\frac{\gamma_m ^2}{r_m^2} \frac{d^2G}{d v^2}.
$$
So the three relationships of $v$ with $s$ that we have take the form
$$
v = \frac{\gamma_m}{r_m} s, \;
\frac{dG}{dv} = \frac{r_m}{\gamma_m} \frac{dG}{ds}, \; \text{and }
\frac{d^2 G}{dv^2} = \frac{r_m^2}{\gamma_m^2} \frac{d^2G}{dv^2}.
$$
Substituting these definitions results in
$$
\gamma_m \left( \frac{\gamma_m}{r_m} s \right) 
\frac{r_m^2}{\gamma_m^2} \frac{d^2 G}{ds^2}
+ \left[ k^{(p)}_{\text{off}} + k^{(p)}_{\text{on}} 
- r_m \left( \frac{\gamma_m}{r_m} s \right) \right] 
\frac{r_m}{\gamma_m} \frac{dG}{ds}
- \frac{r_m k^{(p)}_{\text{on}}}{\gamma_m} G(s) = 0.
$$
Upon simplifying terms, we find an equation that is now in the form of Eq.
$\ref{eq:kummer}$
$$
s \frac{d^2 G}{ds^2}
+ \left(\frac{k^{(p)}_{\text{off}} + k^{(p)}_{\text{off}}}{\gamma_m} - s \right)
\frac{dG}{ds}
- \frac{k^{(p)}_{\text{on}}}{\gamma_m} G(s) = 0.
\label{eq:gen_kummer}
$$

Having put this in the form of the Kummer Eq., we can use Eq.
$\ref{eq:kummer_sol}$ to write $G(s)$ as
$$
\begin{aligned}
G(s) &= A {}_1F_1 
\left(
    \frac{k^{(p)}_{\text{on}}}{\gamma_m}, 
    \frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}{\gamma_m}, 
    s
\right) \\
&+ B s^{1 - \frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}{\gamma_m}}
{}_1 F_1
\left(
    \frac{k^{(p)}_{\text{on}}}{\gamma_m} + 1 -
    \frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}{\gamma_m},
    2 - 
    \frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}{\gamma_m},
    s
\right).
\end{aligned}
$$
We can write down this solution in terms of the original variable of the 
generating function. We have that $s = r_m/\gamma_m v$, and $v = z - 1$. With
this we then write
$$
\begin{aligned}
G(z) &= A {}_1F_1 
\left(
    \frac{k^{(p)}_{\text{on}}}{\gamma_m}, 
    \frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}{\gamma_m}, 
    \frac{r_m}{\gamma_m}(z - 1)
\right) \\
&+ B \left[\frac{r_m}{\gamma_m}(z - 1)\right]
^{1 - \frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}{\gamma_m}}
{}_1 F_1
\left(
    \frac{k^{(p)}_{\text{on}}}{\gamma_m} + 1 -
    \frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}{\gamma_m},
    2 - 
    \frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}{\gamma_m},
    \frac{r_m}{\gamma_m}(z - 1)
\right).
\end{aligned}
\label{eq:gen_sol}
$$

#### Finding the coefficients for the solution

We can now use the normalization condition for the generating function, this is,
$$
G(1) = \sum_{m=0}^\infty 1^m P(m) = 1.
$$
Evaluating $z=1$ in Eq. $\ref{eq:gen_sol}$ results in
$$
G(1) = A {}_1F_1 
\left(
    \frac{k^{(p)}_{\text{on}}}{\gamma_m}, 
    \frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}{\gamma_m}, 
    0
\right).
\label{eq:gen_sol_1}
$$
Let's look at the hypergeometric function evaluated of the form ${}_1F_1(a, b ,
0)$. This takes the form
$$
{}_1 F_{1}(a, b, 0)=\sum_{m=0}^{\infty} \frac{a^{(m)} 0^{n}}{b^{(m)} m!} 
$$
All of the terms but one ($n = 0$) are zero. The first term involving $0^0$ is
undefined. Taking the limit as $z \rightarrow 0$ from the positive side, we find
$$
{}_1F_{1}(a, b, 0) = 
\lim _{z \rightarrow 0^{+}} {}_1 F_{1}(a, b, z) = 
\lim _{z \rightarrow 0^{+}} z^{0} = 1.
$$
Using this property in Eq. $\ref{eq:gen_sol_1}$ tells us that $A = 1$.

We do not have another constraint for $B$. Nevertheless, recall that Eq.
$\ref{eq:first_mom_gen}$ tells us how to compute the first moment of the
distribution from the generating function. For this, we need to compute the
derivative of the confluent hypergeometric function. Let us derive this
identity. Rather than computing the derivative directly, we will compute
$$
z \frac{d}{dz}{}_1F_1 = 
z \frac{d}{dz} 
\left[ \sum_{m=0}^{\infty} \frac{a^{(m)} z^{m}}{b^{(m)} m!}\right].
$$
Taking the derivative inside the sum gives
$$
z \frac{d}{dz}{}_1F_1 = 
z 
\left[ 
    \sum_{m=0}^{\infty} \frac{a^{(m)} }{b^{(m)} m!}
    \frac{d}{dz} z^{m}
\right] 
=
\left[ 
    \sum_{m=0}^{\infty} \frac{a^{(m)} }{b^{(m)}}
    \frac{m z^m}{m!}
\right] .
$$
Simplifying the term $m/m!$ gives
$$
z \frac{d}{dz}{}_1F_1 = 
\left[ 
    \sum_{m=0}^{\infty} \frac{a^{(m)} }{b^{(m)}}
    \frac{z^m}{(m-1)!}
\right] .
\label{eq:confluent_1}
$$
Note that the rising factorials can be rewritten as
$$
\begin{aligned}
    a^{(m)} &=a(a+1)(a+2) \cdots(a+m-1) \\
    &=a \cdot(a+1)[(a+1)+1][(a+1)+2] \cdots[(a+1)+m - 2] \\
    &=a \cdot(a+1)^{(m-1)} .
\end{aligned}
$$
Therefore we can rewrite Eq. $\ref{eq:confluent_1}$ as
$$
\begin{aligned}
\sum_{m=0}^{\infty} \frac{a^{(m)}}{b^{(m)}} \frac{z^{m}}{(m-1) !} &=
\sum_{m=0}^{\infty} \frac{a \cdot(a+1)^{(m-1)}}{b \cdot(b+1)^{(m-1)}} 
\frac{z \cdot z^{(m-1)}}{(m-1) !} \\
&=\frac{a z}{b} 
\sum_{m=0}^{\infty} \frac{(a+1)^{(m-1)}}{(b+1)^{(m-1)}} \frac{z^{m-1}}{(m-1) !}
\end{aligned}
$$
If we define $m' = m - 1$ we have
$$
\frac{a z}{b} \sum_{m=0}^{\infty} \frac{(a+1)^{(m-1)}}{(b+1)^{(m-1)}} 
\frac{z^{m-1}}{(m-1) !} = 
\frac{a z}{b} \sum_{m'=-1}^{\infty}
\frac{(a+1)^{m'}}{(b+1)^{m'}} 
\frac{z^{m'}}{m'!}
$$
The term on the left is almost of the form of the confluent hypergeometric
function again. The only difference is that the sum starts at $m' = -1$. This
first term of the sum would then involve a term of the form $1 / (-1)!$. But
what does this even mean? To find this out, we can generalize the factorial
function using the Gamma function such that
$$
(x - 1)! = \Gamma(x).
$$
The Gamma function diverges as $x \rightarrow 0$, therefore $1/\Gamma(x)
\rightarrow 0$ as $x \rightarrow 0$. This means that the first term of the sum
is zero, so we can begin the sum at $m' = 0$, recovering a confluent
hypergeometric function. With this, we find that
$$
z \frac{d}{d z} {}_1F_1 = 
\frac{a z}{b} \sum_{m=0}^{\infty} \frac{(a+1)^{m}}{(b+1)^{m}} 
\frac{z^{m}}{m !} = 
\frac{a}{b} z_{1} F_{1}(a+1, b+1, z),
$$
therefore
$$
\frac{d}{dz}{}_1F_1 = \frac{a}{b} {}_1F_1(a + 1, b + 1, z).
\label{eq:gen_deriv}
$$

After this small but necessary detour, we can come back to computing the first
moment of our distribution from the generating function. to evaluate Eq.
$\ref{eq:first_mom_gen}$ on Eq. $\ref{eq:gen_sol}$ we first compute the 
derivative of the generating function. This can be easily evaluated using
the relationship we derived for derivatives of ${}_1F_1$. The only thing to be
aware of is that of the chain rule. In particular for our third entry of the
third entry of the function we have $r_m / \gamma_m (z - 1)$ rather than simply
$z$ as we had in Eq. $\ref{eq:gen_deriv}$. This means that by the chain rule we
have that if we define $u = r_m / \gamma_m (z - 1)$, we have
$$
du = \frac{r_m}{\gamma_m} dz \Rightarrow
\frac{dG}{dz} = 
\frac{dG}{du} \frac{du}{dz} = 
\frac{dG}{du} \frac{r_m}{\gamma_m}.
$$
So there is an extra factor of $r_m / \gamma_m$ that will come along when we
compute the derivative of our generating functions. Computing the derivative
of Eq. $\ref{eq:gen_sol}$ results in
$$
\small
\begin{aligned}
\frac{dG}{dz} &= 
\frac{k^{(p)}_{\text{on}}}{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}
\frac{r_m}{\gamma_m}
{}_1F_1 
\left(
    \frac{k^{(p)}_{\text{on}}}{\gamma_m} + 1, 
    \frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}{\gamma_m} + 1, 
    \frac{r_m}{\gamma_m}(z - 1)
\right) \\
& + B 
\left( 1 - \frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}{\gamma_m} \right)
\left[\frac{r_m}{\gamma_m}(z - 1)\right]
^{\frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}{\gamma_m}}
{}_1 F_1
\left(
    \frac{k^{(p)}_{\text{on}}}{\gamma_m} + 1 -
    \frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}{\gamma_m},
    2 - 
    \frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}{\gamma_m},
    \frac{r_m}{\gamma_m}(z - 1)
\right) \\
&+ B \left[\frac{r_m}{\gamma_m}(z - 1)\right]
^{1 - \frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}{\gamma_m}}
\left( 
    \frac{k^{(p)}_{\text{on}} + \gamma_m}
    {k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}} + \gamma_m}
\right)
\frac{r_m}{\gamma_m}
{}_1 F_1
\left(
    \frac{k^{(p)}_{\text{on}}}{\gamma_m} + 2 -
    \frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}{\gamma_m},
    1 - 
    \frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}{\gamma_m},
    \frac{r_m}{\gamma_m}(z - 1)
\right).
\end{aligned}
$$
This rather convoluted result is enormously simplified upon evaluating the
derivative at $z = 1$ (See Eq. $\ref{eq:first_mom_gen}$). This results in
$$
\left. \frac{dG}{dz} \right\vert_{z = 1} = 
\frac{k^{(p)}_{\text{on}}}{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}
\frac{r_m}{\gamma_m}
{}_1F_1 
\left(
    \frac{k^{(p)}_{\text{on}}}{\gamma_m} + 1, 
    \frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}{\gamma_m} + 1, 
    0
\right) = 
\frac{r_m}{\gamma_m}
\frac{k^{(p)}_{\text{on}}}{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}},
$$
which is exactly the mean mRNA copy number we derived before. Since $B$ does not
contribute to the mean, we can safely assume that $B = 0$. This means that
the final result for the generating function takes the much more compact form
$$
G(z) = 
{}_1F_1 
\left(
    \frac{k^{(p)}_{\text{on}}}{\gamma_m}, 
    \frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}{\gamma_m}, 
    \frac{r_m}{\gamma_m}(z - 1)
\right).
\label{eq:gen_final}
$$

#### Extracting the steady-state mRNA distribution

Let us quickly recapitulate where we are. We started with a system of infinite
many ordinary differential equations, one for each promoter state and mRNA copy
number that defined the master equation for our two-state promoter. We then used
the generating function to transform this system into a single partial
differential equation. The resulting differential equation for the generating
function took the form of the so-called Kummer differential equation, which has
as a solution the confluent hypergeometric function and the Tricomi function.
After imposing the normalization condition on the generating function, we found
that the confluent hypergeometric function's coefficient was $A=1$. We then used
the fact that the mean mRNA copy number $\langle m \rangle$ exists to show that
the Tricomi function's coefficient is $B=0$. All that effort lead us to Eq.
$\ref{eq:gen_final}$, the generating function for the two-state promoter mRNA
steady-state distribution. All we have left is trying to beat Eq.
$\ref{eq:gen_final}$ into the form of a standard generating function to extract
the probability distribution from it.

Let us begin this task by writing down Eq. $\ref{eq:gen_final}$ with the full
definition of the confluent hypergeometric function. This gives us
$$
G(z) = \sum_{m=0}^\infty 
\frac{
    \left(\frac{k^{(p)}_{\text{on}}}{\gamma_m}\right)^{(m)}
}{
    \left(\frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}
    {\gamma_m}\right)^{(m)}
}
\frac{
    \left[\frac{r_m}{\gamma_m} (z-1) \right]^m
}{
    m!
}
$$
Let us now split apart the term $(z-1)$, obtaining
$$
G(z) = \sum_{m=0}^\infty 
\frac{
    \left(\frac{k^{(p)}_{\text{on}}}{\gamma_m}\right)^{(m)}
}{
    \left(\frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}
    {\gamma_m}\right)^{(m)}
}
\frac{
    \left(\frac{r_m}{\gamma_m} \right)^m
}{
    m!
}
(z - 1)^m.
$$
We now rewrite this last term $(z-2)^m$ using the binomial expansion. This
results in
$$
G(z) = \sum_{m=0}^\infty 
\frac{
    \left(\frac{k^{(p)}_{\text{on}}}{\gamma_m}\right)^{(m)}
}{
    \left(\frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}
    {\gamma_m}\right)^{(m)}
}
\frac{
    \left(\frac{r_m}{\gamma_m} \right)^m
}{
    m!
}
\left[ 
    \sum_{n=0}^m {m \choose n} z^n (-1)^{m -  n}.
\right]
$$
We can take out the sum over the index $n$ to the front, obtaining
$$
G(z) = \sum_{m=0}^\infty \sum_{n=0}^n
\frac{
    \left(\frac{k^{(p)}_{\text{on}}}{\gamma_m}\right)^{(m)}
}{
    \left(\frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}
    {\gamma_m}\right)^{(m)}
}
\frac{
    \left(\frac{r_m}{\gamma_m} \right)^m
}{
    m!
}
\left[ 
    {m \choose n} z^n (-1)^{m -  n}
\right].
$$
To make further progress, we must reindex the sum. The trick is to reverse the
default order of the sums as
$$
\sum_{m=0}^{\infty} \sum_{n=0}^{m} = \sum_{n=0}^{\infty} \sum_{m=n}^{\infty}.
$$
To see the logic of the sum, we point the reader to [@Fig:ch5_fig37]. The key is
to notice that the double sum $\sum_{m=0}^\infty \sum_{n=0}^m$ is adding all
possible pairs $(m, n)$ in the lower triangle, so we can add the terms
vertically as the original sum indexing suggests, i.e.
$$
\sum_{m=0}^{\infty} \sum_{n=0}^{m} x_{(m, n)}= 
x_{(0, 0)} + x_{(1, 0)} + x_{(1, 1)} + x_{(2, 0)} + x_{(2, 1)} + x_{(2, 2)} + 
\ldots,
$$
where the variable $x$ is just a placeholder to indicate the order in which the
sum is taking place. But we can also add the terms horizontally as
$$
\sum_{n=0}^{\infty} \sum_{m=n}^{\infty} x_{(m, n)} =
x_{(0, 0)} + x_{(1, 0)} + x_{(2, 0)} + \ldots + x_{(1,1)} + x_{(2, 1)} + \ldots,
$$
which still adds all of the lower triangle terms.

![**Reindexing double sum.** Schematic for reindexing the sum $\sum_{m=0}^\infty
\sum_{n=0}^m$. Blue circles depict the 2D grid of nonnegative integers
restricted to the lower triangular part of the $m, n$ plane. The trick is that
this double sum runs over all $(m, n)$ pairs with $n\le m$. Summing $m$ first
instead of $n$ requires determining the boundary: the upper boundary of the
$n$-first double sum becomes the lower boundary of the $m$-first double
sum.](ch5_fig37){#fig:ch5_fig37 short-caption="Reindexing double sum"}

Rewriting the sum in this way results in
$$
G(z) = \sum_{n=0}^\infty \sum_{m=n}^\infty
\frac{
    \left(\frac{k^{(p)}_{\text{on}}}{\gamma_m}\right)^{(m)}
}{
    \left(\frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}
    {\gamma_m}\right)^{(m)}
}
\frac{
    \left(\frac{r_m}{\gamma_m} \right)^m
}{
    m!
}
\left[ 
    {m \choose n} z^n (-1)^{m -  n}
\right].
$$
This allows us to separate the variable $z^n$ from the rest of the equation,
leaving the standard format generating function to read the probability
distribution $P(m)$. This looks as
$$
G(z) = \sum_{n=0}^\infty z^n 
\left[
\sum_{m=n}^\infty
\frac{
    \left(\frac{k^{(p)}_{\text{on}}}{\gamma_m}\right)^{(m)}
}{
    \left(\frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}
    {\gamma_m}\right)^{(m)}
}
\frac{
    \left(\frac{r_m}{\gamma_m} \right)^m
}{
    m!
}
    {m \choose n} (-1)^{m -  n}
\right].
$$
Given the "dummy" nature of $z$, it does not matter what the sum variable name 
is. We can simply rename $m = n$ and $n = m$ and conclude that our distribution
takes the form
$$
P(m) = 
\sum_{n=m}^\infty
\frac{
    \left(\frac{k^{(p)}_{\text{on}}}{\gamma_m}\right)^{(n)}
}{
    \left(\frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}
    {\gamma_m}\right)^{(n)}
}
\frac{
    \left(\frac{r_m}{\gamma_m} \right)^n
}{
    n!
}
    \frac{n!}{m! (n - m)!} (-1)^{n -  m}.
\label{eq:prob_ss_1}
$$

We can simplify Eq. $\ref{eq:prob_ss_1}$ further. First we split the term 
$(-1)^{n-m} = (-1)^{-m} (-1)^{n}$. Furthermore we absorbed the $(-1)^{n}$ term
on the $(r_m / \gamma_m)^n$ term. We also cancel the obvious $n!/n!$ term, 
obtaining
$$
P(m) = \sum_{n = m}^\infty
\frac{(-1)^{-m}}{m!}
\frac{
    \left( \frac{k^{(p)}_{\text{on}}}{\gamma_m}\right)^{(n)}
}{
    \left( \frac{k^{(p)}_{\text{on}}+ k^{(p)}_{\text{off}}}
    {\gamma_m}\right)^{(n)}
}
\frac{\left( - \frac{r_m}{\gamma_m}\right)^n}{(n - m)!}.
\label{eq:prob_ss_2}
$$
We recognize in Eq. $\ref{eq:prob_ss_2}$ that we have almost all the terms for a
confluent hypergeometric function ${}_1F_1$. The problem is that the sum starts
at $n = m$ rather than $n = 0$. Since the upper limit of the sum is $\infty$, we
can simply define $u = n - m \Rightarrow n = m + u$. We can then use the 
following property of raising factorials
$$
\begin{split}
a^{(n)} &=a(a+1)(a+2) \cdots(a+n-1), \\
&=a(a+1)(a+2) \cdots(a+(u+m)-1), \\
&=a(a+1) \cdots(a+m-1)(a+m)(a+m+1) \cdots(a+m+u-1), \\
&=a^{(m)}(a+m)^{(u)}.
\end{split}
$$
Making these substitutions results in
$$
P(m) = \sum_{u=0}^{\infty} 
\frac{(-1)^{-m}}{m !} 
\frac{
    \left(\frac{k^{(p)}_{\text{on}}}{\gamma_m} \right)^{(m)}
    \left(\frac{k^{(p)}_{\text{on}}}{\gamma_m} + m \right)^{(u)}
    \left(-\frac{r_m}{\gamma_m}\right)^{u}
    \left(-\frac{r_m}{\gamma_m}\right)^{m}
}{
    \left(
        \frac{
            k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}
        }{
            \gamma_m
        }
    \right)^{(m)}
    \left(
        \frac{
            k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}
            }{
                \gamma_m
            } + m
    \right)^{(n)}}
    \frac{1}{u!}.
$$
Taking out of the sum the terms that do not depend on $u$ gives
$$
P(m) = 
\frac{(-1)^{-m}}{m!}
\frac{
    \left(
        \frac{
            k^{(p)}_{\text{on}}
        }{
            \gamma_m
        }
    \right)^{(m)}
}{
    \left(
        \frac{
            k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}
        }{
            \gamma_m
        }
    \right)^{(m)}
}
\left(- \frac{r_m}{\gamma_m}\right)^m
\left[
    \sum_{u=0}^{\infty}
    \frac{
        \left(
            \frac{
                k^{(p)}_{\text{on}} 
            }{
                \gamma_m
            }
            + m
        \right)^{(u)}
    }{
        \left(
            \frac{
                k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}} 
            }{
                \gamma_m
            }
            + m
        \right)^{(u)}
    }
    \frac{
        \left(- \frac{r_m}{\gamma_m}\right)^u
    }{u!}
\right].
$$
We recognize the term in the square brackets to be the necessary components for
a confluent hypergeometric function. We can therefore write the mRNA
steady-state distribution as
$$
P(m) = 
\frac{1}{m!}
\frac{
    \left(
        \frac{
            k^{(p)}_{\text{on}}
        }{
            \gamma_m
        }
    \right)^{(m)}
}{
    \left(
        \frac{
            k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}
        }{
            \gamma_m
        }
    \right)^{(m)}
}
\left(\frac{r_m}{\gamma_m}\right)^m
{}_1F_1 
\left(
    \frac{
            k^{(p)}_{\text{on}} 
        }{
            \gamma_m
        }
    + m,
    \frac{
            k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}} 
        }{
            \gamma_m
        }
    + m,
    - \frac{r_m}{\gamma_m}
\right).
$$
For the last ingredient, we remove the rising factorials using the identity
$$
\begin{split}
a^{(m)} &=(a)(a+1)(a+2) \cdots(a+m-1), \\
&=\frac{(a+m-1) \cdots(a)(a-1) \cdots (1)}{(a+1) \cdots(1)}, \\
&=\frac{(a+m-1) !}{(a-1) !}.
\end{split}
$$
This allows us to write
$$
\begin{split}
P(m) &= 
\frac{1}{m!}
\frac{
    \left(
        \frac{k^{(p)}_{\text{on}}}{\gamma_m} + m - 1
    \right) !
}{
    \left(
        \frac{k^{(p)}_{\text{on}}}{\gamma_m} - 1
    \right) !
}
\frac{
    \left(
        \frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}{\gamma_m} - 1
    \right) !
}{
    \left(
        \frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}{\gamma_m} + m - 1
    \right) !
}
\left( \frac{r_m}{\gamma_m} \right)^m \\
&\times {}_1F_1 
\left(
    \frac{
            k^{(p)}_{\text{on}} 
        }{
            \gamma_m
        }
    + m,
    \frac{
            k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}} 
        }{
            \gamma_m
        }
    + m,
    - \frac{r_m}{\gamma_m}
\right).
\end{split}
$$
Or in terms of Gamma functions, we obtain the final form of the steady-state
mRNA distribution
$$
\begin{aligned}
P(m) &= 
\frac{1}{\Gamma(m + 1)}
\frac{
    \Gamma
    \left(
        \frac{k^{(p)}_{\text{on}}}{\gamma_m} + m
    \right)
}{
    \Gamma
    \left(
        \frac{k^{(p)}_{\text{on}}}{\gamma_m}
    \right)
}
\frac{
    \Gamma
    \left(
        \frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}{\gamma_m}
    \right)
}{
    \Gamma
    \left(
        \frac{k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}}}{\gamma_m} + m 
    \right)
}
\left( \frac{r_m}{\gamma_m} \right)^m \\
&\times {}_1F_1 
\left(
    \frac{
            k^{(p)}_{\text{on}} 
        }{
            \gamma_m
        }
    + m,
    \frac{
            k^{(p)}_{\text{on}} + k^{(p)}_{\text{off}} 
        }{
            \gamma_m
        }
    + m,
    - \frac{r_m}{\gamma_m}
\right),
\end{aligned}
$$
The equation used to fit the kinetic parameters for the unregulated promoter.
