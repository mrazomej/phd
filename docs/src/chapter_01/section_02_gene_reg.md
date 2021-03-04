## Gene regulation as a Physics 101 problem

As organisms navigate the challenges presented by the environment, they must 
constantly fight against the will of the second law of thermodynamics to bring
them back to an equilibrium state. To face those challenges, cells are equipped
with a toolkit of genes written in the language of A, T, C, and G of the genome.
We can think of a typical bacteria genome with $\approx 5\times 10^3$ genes as 
the instruction manual to produce a repertoire of tools that allow cells to
thrive under different circumstances that they face throughout their lives. But
not all challenges are the same, and therefore not all tools are needed 
simultaneously. Therefore all organisms are faced with the challenge of 
orchestrating the expression of the correct subset molecular tools at their
disposal when trying to survive on different environments. From cells in the
fly embryo expressing different genes that will define their identity on the
final body plan of the animal, to a simple bacteria expressing the correct
enzymes to process the available nutrients in the environment.

Our understanding of how organisms regulate the expression of their genes is
still not as thorough as one might expect given the amount of effort that has
gone into this question. Take for example *E. coli*--arguably the most well
characterized model organism--for which we know the regulatory scheme of less
than 1/3 of its genes [@Ireland2020]. For more complex organisms such as
*Drosophila*, *C. elegans*, or even humans we are even more hopeless on getting
a holistic view of the regulatory landscape. Nevertheless we would not be doing
justice to the great advances in the field if we were to pretend we are
completely ignorant about how gene regulation takes place in bacteria. There is
a rich mechanistic understanding of how the transcriptional machinery takes the
information contained in the DNA and transcribes it into RNA [@Browning2004].
The relative simplicity of the process has inspired generations of biophysicists
to try to write down minimal models that can describe and predict features of
the process of gene regulation [@Ackers1982;@Bintu2005;@Kuhlman2007].

These modeling efforts come into two main flavors: equilibrium statistical 
mechanical models, and kinetic models. In the following sections we will 
introduce the necessary background for both approaches relevant to the rest of
the thesis.

### Minimal model of gene expression

Let us begin our introduction to gene expression modeling with the simplest
example. As shown in [@Fig:ch1_fig01](A) we imagine a gene promoter (the region
of the gene from which transcriptional regulation takes place) produces mRNA at
a constant rate $r_m$. Each individual mRNA can stochastically decay with a rate
$\gamma_m$. Our interest is to understand how the mRNA count $m$ changes over
time given these two competing processes. For that let us write the mRNA count
at time $m(t + \Delta t)$, where $t$ is the time--which we are thinking of as
being "right now"--and $\Delta t$ is a tiny time step into the future. The mRNA
count can then be predicted by computing
$$
m(t + \Delta t) = m(t) + r_m \Delta t - (\gamma_m \Delta t) m(t),
$$
where we can think of $r_m \Delta t$ as the probability of observing a single
mRNA being produced in the time interval $[t, t + \Delta t]$ ($\Delta t$ is so
small that we neglect the possibility of seing multiple mRNAs being produced),
and $\gamma_m \Delta t$ the probability of seeing a single mRNA being degraded.
But since each mRNA has the same probability of being degraded, the total number
of mRNAs that we would see decay in this time window would be the probability
per mRNA times the total number of mRNAs. If we send the term $m(t)$ to the left
hand side of the equation and divide both sides by $\Delta t$, we obtain
$$
\frac{m(t + \Delta t) - m(t)}{\Delta t} = r_m - \gamma_m m(t).
$$
Upon taking the limit when $\Delta t \rightarrow 0$ we see that the left hand
side is the definition of the derivative of the mRNA count with respect to time.
We then obtain an ordinary differential equation of the form
$$
\frac{dm}{dt} = r_m - \gamma_m m(t).
\label{eq:dm_dt}
$$
Before even attempting to solve $\ref{eq:dm_dt}$ we can perform a qualitative
analysis of the dynamics [@Strogatz2018]. In particular it is useful to plot the
contribution to the derivative $dm/dt$ for each of the components (production
and degradation) as a function of $m$. This is shown in [@Fig:ch1_fig01](B)
where the blue horizontal line $r_m$ shows the production rate--which does not
depend on $m$, and the red line shows the degradation term $m\gamma_m$ which
scales linearly with $m$. Notice that for the degradation term we are not
including the negative sign, i.e., we are not plotting $-m \gamma_m$. The point
$m_{ss}$ where both lines intersect represent the point where the production
matches the degradation. For all values before $m_{ss}$ the production term is
larger than the degradation, which means that for any value $m < m_{ss}$ the
derivative is positive ($dm/dt > 0$), so over time the system will produce more
mRNA. The opposite is true for values all values after $m_{ss}$ where the
degradation term is larger than the production term, implying that $dm/dt < 0$.
This means that for $m > m_{ss}$ the system will degrade mRNA. This opposite
trends point at the idea that $m_{ss}$ must be what is called a stable fix point
of the dynamical system. This can schematically be seen at the bottom of
[@Fig:ch1_fig01](B) where the size of the arrow heads indicates the trend of the
system to move either left or right in $m$. Since all arrows point at the
special value $m_{ss}$ we can say that any small perturbation of the system will
be dissipated as the system relaxes back to $m_{ss}$.

This qualitative statement can be confirmed by solving Eq. $\ref{eq:dm_dt}$. If
we define the initial condition $m(t=0) = m_o$ by separation of variables we 
will obtain a solution of the form
$$
m(t) = m_o e^{-\gamma_m t} + \frac{r_m}{\gamma_m} (1 - e^{-\gamma_m t}).
$$
In the limit when $t \rightarrow \infty$ we can see that the steady state 
solution is given by
$$
m_{ss} = \frac{r_m}{\gamma_m}.
$$
[@Fig:ch1_fig01](C) shows the time evolution of $m$ for different initial values
$m_o$. We can see that indeed regardless of the initial mRNA count the system
relaxes exponentially fast to $m_{ss} = r_m / \gamma_m$.

![**Minimal model of gene expression.** (A) Schematic of the kinetics governing
gene expression. mRNA is produced at a constant rate $r_m$ independent of the
current mRNA copy number. Degradation of each mRNA occurs at a rate $\gamma_m$.
(B) Example of the qualitative analysis of the mRNA dynamics via a 1D
phase-portrait. The differential equation governing the dynamics contains two
terms: a constant production rate given by $r_m$, and a degradation rate
$\gamma_m m$ that depends on the current mRNA count. The main plot shows each of
the components in the $m$ vs $dm/dt$ plot. Since $r_m$ does not depend on the
current number of mRNA it gives a straight production rate as a function of $m$.
The total degradation rate depends linearly with the mRNA copy number, giving a
line with slope $\gamma_m$. When the two components are equal (bot lines
crossing), we obtain the steady-state mRNA value $m_{ss}$. The bottom line shows
a qualitative schematic of the flow of the system towards this steady state. The
further $m$ is from $m_{ss}$, the faster it moves towards this point as
schematized by the size of the arrows. (C) Example of mRNA dynamics for
different initial conditions. Over time all curves converge to the steady-state
mRNA value $m_{ss}=r_m/\gamma_m$. For this plot $\gamma_m = 1$ and $r_m/\gamma_m
= 10$. The [Python code
(`ch1_fig01C.py`)](https://github.com/mrazomej/phd/blob/master/src/chapter_01/code/ch1_fig01C.py)
used to generate part (C) of this figure can be found on the thesis [GitHub
repository](https://github.com/mrazomej/phd).](ch1_fig01){#fig:ch1_fig01
short-caption="Minimal model of gene expression"}

So far our model assumes a simple constant transcription rate $r_m$; let us
expand this term a little further in order to include regulation by a
transcriptional repressor further down the road. We knot that for a
transcriptional event to take place the RNA polymerase (RNAP) must bind to the
promoter region and undergo a series of irreversible steps such as opening the
double helix to initiate the copying of the DNA sequence into mRNA
[@Browning2004]. If we assume these irreversible steps take place in a much 
longer time-scale compared to the initial binding and unbinding of the RNAP on
the promoter, we can separate these two events in different terms. In particular
we can write that the mRNA production happens with a rate
$$
\text{mRNA production} = r_m \cdot p_{\text{bound}},
$$
where we split the original production term into two steps: $p_{\text{bound}}$,
the probability of finding an RNAP bound to the promoter, and $r_m$ which
captures all of the downstream irreversible steps that take place once the RNAP
is engaged in a transcriptional event. A way to think about it--relevant to what
I am doing right now as I type my thesis--is to think that the speed at which I
type this document has to do with two things: The probability of me being
actively working on these notes times the rate at which I type these notes once
I engage in the activity. The reason this splitting makes sense is because we
can include the effect of the regulation by a transcriptional repressor as a
reduction of the time (the probability) that the RNAP can be bound to the
promoter by effectively taking away probability mass. Futhermore, since we are
assuming that the binding and unbinding of the RNAP happen at a timescale much
faster than the downstream events, we can assume this binding reaction is in
quasi-equilibrium, for which we can use the powerful theoretical framework of
statistical mechanics. Let us now delve into the basics of this physical theory.

### The unreasonable effectiveness of unrealistic simplifications

On the preface of the textbook *Molecular Driving Forces* Dill and Bromberg
introduce the idea of Statistical Mechanics as *the unreasonable effectiveness
of unrealistic simplifications* [@Dill2010]. Although one could make the case
that all of physics follows this description, it is certainly evident that
statistical mechanics is a vivid example of how simple ideas can have profound
consequences. Statistical mechanics can be simply defined as the theory that,
upon assuming the atomic nature of matter, explains the phenomenology that
classical thermodynamics established [@Dill2010] from the interactions of the
microscopic components of a system. As with any other physical theory,
statistical mechanics is built from a set of *empirical* facts that define
"axioms" that we take to be truth. In other words, as Feynman famously described
to us: if we want to come up with a new law of nature there is a simple recipe
that we must follow:

1. We guess the law. Literally. The deepest understanding of physical reality we
   have up to now come from educated guesses made after a careful observation of
   nature.

2. We compute the consequences of such guess. That is why mathematical theories
   allow us to sharpen our statements about how we think nature works.

3. We compare with experiments/observations. The scientific revolution came
   about when, after the dark ages, we finally learned it was okay to say "we
   don't know."

In such simple statement, Feynman tells us, lies the key to science
[@Feynman1965]. For our purpose of understanding the basis of statistical
mechanics we will argue that the main law upon which the field is founded is
given by Boltzmann's law
$$
\frac{P(E_1)}{P(E_2)} = \frac{e^{-E_1 / k_BT}}{e^{-E_2 / k_BT}}.
\label{eq:boltzmann_law}
$$
Let us unpack this equation. The main idea behind statistical mechanics is that
macroscopic observables (temperature, pressure, specific heat in classic
examples) are emergent properties dictated by the dynamics of the microscopic
components of the system. What Boltzmann's law tells us is that the relative
probability of a system in thermal equilibrium to be found in a particular
microstate with energy $E_1$ compared to being in a microstate with energy $E_2$
is given by an exponential function of minus the energy of such microstate
divided by $k_BT$, the thermal energy. To give concrete examples of what a
microstate can look like, [@Fig:ch1_fig02](A) shows three molecular systems
relevant for biology. On the first example we have the classic ligand-receptor
binding problem; here we imagine that a solution can be discretized into a
series of small boxes. In each of these boxes one and only one ligand molecule
can fit in. In principle we can list all possible spatial arrangements of
ligands. We could then calculate the relative likelihood of finding the system
in any of the configurations as long as we can assign an energy value to each of
them. The second example focuses on ligand-gated ion channels. In this
particular system we care about the state of the ion channel itself--either open
or close--and the binding configuration of the ligands. If the channel responds
to the concentration of the ligand by changing its probability of gating, that
is something we can calculate using equilibrium statistical mechanics. Finally,
the third example shows different configurations of a small patch of cell
membrane. All deformations of a membrane have energetic costs associated with
them. So listing all possible membrane configurations we can calculate what is
the most likely shape of a membrane given the forces and stresses acting on it. 

The macroscopic states that we get to observe can then be thought of as a
coarse-graining of many microstates into a single macrostate. For example, in
the case of the ligand-receptor binding, we rarely would care about the specific
position of all the ligand molecules in the solution. What we would be
interested in is whether or not the ligand is bound to the receptor. We can
therefore define as our "macrostate" the particular configuration of the
receptor as schematically shown in [@Fig:ch1_fig02](B).

![**Boltzmann's law and the definition of a micro and macrostate.** (A) Top
panel: ligand-receptor binding microstates. Middle panel: ligand-gated ion
channel microstates. Bottom panel: membrane patch deformations. (B) Schematic of
the definition of a "macrostate." In the ligand-receptor binding problem we
ignore the spatial configuration of all ligand molecules, and focus on the
binding state of the receptor.](ch1_fig02){#fig:ch1_fig02
short-caption="Boltzmann's law and the definition of a micro and macrostate"}

If we want to know the likelihood of finding a particular system in any specific
configuration Boltzmann's law (Eq. $\ref{eq:boltzmann_law}$) is then telling us
a protocol we must follow: 

1. Enumerate all possible microstates in which the system can be found.

2. Compute the energy of each of these microstates.

3. Define the "macrostate" we care about by grouping all microstaes that belong
   to the same macrostate.

4. Compute the Boltzmann factor by exponentiating minus the energy divided by
   the thermal energy.

To see this protocol in action let us apply it to the calculation of
$p_{\text{bound}}$, the probability of finding an RNAP bound to the promoter. We
will go through each of the steps on the protocol and build up the "unrealistic
simplifications" that will allow us to make this calculation.

**1. Enumerate possible microstates.** We begin by making a drastic
coarse-graining of the bacterial genome. For us a genome is simply made out of
boxes where the RNAP can bind. We imagine that there is a single site where RNAP
can bind specifically--the promoter of interest. There are also $N_{NS} \approx
5\times 10^6$ non-specific binding sites, one per basepair (bp) in the genome.
We ignore the fact that the RNAP footprint when it binds to the genome is
roughly 30 bp. This assumption is valid if the number of available RNAP
molecules is much smaller than the number of non-specific binding sites since it
is extremely unlikely that by pure chance two RNAPs would fall next to each
other. We also ignore the possibility of RNAP not being bound to the genome.
Experimental evidence with mini-cells in which a mutant *E. coli* that sheds
vesicles without segregating a chromosome support this assumption [@Bintu2005].
The exercise then consists of choosing at random one box for each of the $P$
polymerase available to bind. [@Fig:ch1_fig03] shows in the first column two
possible configurations of our coarse-grained genome.

**2. Compute the energy for each microstate.** Let us analyze first the case
where all $P$ RNAP molecules are bound non-specifically to the genome. For
simplicity we assume that RNAP binds to all $N_{NS}$ non-specific binding sites
with the same affinity. We assign this energy to be $\varepsilon_P^{(NS)}$. This
assumption could be relaxed as we explored in [@Phillips2019], but for now we
don't have to worry about this complication. For a statistical mechanician the
assignment of binding energies does not come from some quantum first-principled
calculation or anything similar. We simply label the interaction of the RNAP and
the rest of the genome with a single value, $\varepsilon_P^{(NS)}$, that
coarse-grains all of the hydrogen bonds and other effects that go into this
physical process. Since we have $P$ such polymerases bound non specifically, the
energy of any state with a similar configuration is then $P
\varepsilon_P^{(NS)}$ as shown in [@Fig:ch1_fig03] second column, top row.

**3. Define the "macrostate" we care about.** In a sense when we speak about
macrostate, it does not necessarily mean something that we can macroscopically
observe. What it means is that we group togater--a form of coarse-graining--a
bunch of states that we take to be functionally equivalent.

![**Statistical Mechanics protocol for RNAP binding.** On a discretized genome
we follow the statistical mechanics protocol to compute the Boltzmann weight of
each of the relevant microstates. The $P$ available RNAPs are assumed to have
two binding configurations: One specific binding to the promoter of interest
(with energy $\varepsilon_P^{(S)}$) and non-specific to any of the $N_{NS}$
non-specific binding sites (with energy
$\varepsilon_P^{(NS)}$).](ch1_fig03){#fig:ch1_fig03 short-caption="Statistical
Mechanics protocol for RNAP binding"}