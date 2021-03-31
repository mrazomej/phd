## Gene regulation as a Physics 101 problem

As organisms navigate the challenges presented by the environment, they must
constantly fight against the will of the second law of thermodynamics to bring
them back to an equilibrium state. To face such challenges, cells are equipped
with a toolkit of genes written in the language of A, T, C, and G of the genome.
We can think of a typical bacteria genome with $\approx 5\times 10^3$ genes as
the blueprint to produce a repertoire of tools that allow cells to thrive under
myriad circumstances that they face throughout their lives. Given the vast
number of challenges that organisms face, there is constant pressure on every
living system to use the right tools for the right circumstances. Thus all
organisms are faced with the task of orchestrating the expression of the correct
subset of genes at their disposal when trying to survive. From cells in the fly
embryo expressing different genes that will define their identity on the
animal's final body plan to a simple bacteria expressing the correct enzymes to
process the available nutrients in the environment.

Our understanding of how organisms regulate their genes' expression is still not
as thorough as one might expect, given the effort that has gone into this
question. Take, for example, *E. coli*--arguably the most well-characterized
model organism--for which we know the regulatory scheme of less than 1/3 of its
genes [@Ireland2020]. For more complex organisms such as *Drosophila*, *C.
elegans*, or even humans, we are even more hopeless on getting a holistic view
of the regulatory landscape. Nevertheless, we would not be doing justice to the
field's significant advances if we were to pretend we are utterly ignorant about
how gene regulation takes place in bacteria. There is a rich mechanistic
understanding of how the transcriptional machinery takes the information
contained in DNA and transcribes it into RNA [@Browning2004]. The relative
simplicity of the process has inspired generations of biophysicists to try to
write down minimal models that can describe and predict features of the process
of gene regulation [@Ackers1982;@Bintu2005;@Kuhlman2007].

These modeling efforts come into two main flavors: equilibrium statistical
mechanical models and kinetic models. In the following sections, we will
introduce the necessary background for both approaches relevant to the rest of
the thesis.

### Minimal model of gene expression

Let us begin our introduction to gene expression modeling with the simplest
example. As shown in [@Fig:ch1_fig01](A), we imagine a gene promoter (the region
of the gene where transcriptional regulation takes place) produces mRNA at
a constant rate $r_m$. Each mRNA can stochastically decay with a rate
$\gamma_m$. Our interest is to understand how the mRNA count $m$ changes over
time, given these two competing processes. For that, let us write the mRNA count
at time $m(t + \Delta t)$, where $t$ is the time--which we are thinking of as
being "right now"--and $\Delta t$ is a little time step into the future. The
mRNA count can then be predicted by computing
$$
m(t + \Delta t) = m(t) + r_m \Delta t - (\gamma_m \Delta t) m(t),
\label{eq:m_t_Delta_t}
$$
where we can think of $r_m \Delta t$ as the probability of observing a single
mRNA being produced in the time interval $[t, t + \Delta t]$ ($\Delta t$ is so
small that we neglect the possibility of seeing multiple mRNAs being produced),
and $\gamma_m \Delta t$ the probability of seeing a single mRNA being degraded.
But since each mRNA has the same probability of being degraded, the total number
of mRNAs that we would see decay in this time window would be the probability
per mRNA times the total number of mRNAs. This is in contrast with the
production of mRNA, which does not depend on the current number of mRNAs. If we
send the term $m(t)$ to the left-hand side of the equation and divide both sides
by $\Delta t$, we obtain
$$
\frac{m(t + \Delta t) - m(t)}{\Delta t} = r_m - \gamma_m m(t).
$$
Upon taking the limit when $\Delta t \rightarrow 0$, we see that the left-hand
side is the definition of the derivative of the mRNA count with respect to time.
We then obtain an ordinary differential equation of the form
$$
\frac{dm}{dt} = r_m - \gamma_m m(t).
\label{eq:dm_dt}
$$
Before even attempting to solve $\ref{eq:dm_dt}$, we can perform a qualitative
analysis of the dynamics [@Strogatz2018]. It is handy to plot the contribution
of each of the components (production and degradation) to the derivative $dm/dt$
as a function of $m$. This is shown in [@Fig:ch1_fig01](B), where the blue
horizontal line $r_m$ shows the production rate--which does not depend on $m$,
and the red line shows the degradation term $m\gamma_m$ which scales linearly
with $m$. Notice that we do not include the negative sign for the degradation
term, i.e., we are not plotting $-m \gamma_m$. The point $m_{ss}$ where both
lines intersect represents the point where the production matches the
degradation. For all values less than $m_{ss}$ the production term is larger
than the degradation, which means that for any value $m < m_{ss}$ the derivative
is positive ($dm/dt > 0$), so over time the system will accumulate more mRNA.
The opposite is true for all values after $m_{ss}$ where the degradation term is
larger than the production term, implying that $dm/dt < 0$. This means that for
$m > m_{ss}$, the system will tend to lose mRNA. These opposite trends point to
the idea that $m_{ss}$ must be called a stable fixed point of the dynamical
system. This can schematically be seen at the bottom of [@Fig:ch1_fig01](B). The
arrowheads' size indicates the system's trend to move either left or right in
$m$. Since all arrows point at the special value, $m_{ss}$, we can say that any
small perturbation of the system will be dissipated as the system relaxes back
to $m_{ss}$.

This qualitative statement can be confirmed by solving Eq. $\ref{eq:dm_dt}$. If
we define the initial condition $m(t=0) = m_o$ by separation of variables we 
will obtain a solution of the form
$$
m(t) = m_o e^{-\gamma_m t} + \frac{r_m}{\gamma_m} (1 - e^{-\gamma_m t}).
$$
In the limit when $t \rightarrow \infty$ we can see that the steady-state
solution is given by
$$
m_{ss} = \frac{r_m}{\gamma_m}.
$$
[@Fig:ch1_fig01](C) shows the time evolution of $m$ for different initial values
$m_o$. We can see that indeed regardless of the initial mRNA count the system
relaxes exponentially to $m_{ss} = r_m / \gamma_m$.

![**Minimal model of gene expression.** (A) Schematic of the kinetics governing
gene expression. mRNA is produced at a constant rate $r_m$ independent of the
current mRNA copy number. Degradation of each mRNA occurs at a rate $\gamma_m$.
(B) Example of the qualitative analysis of the mRNA dynamics via a 1D
phase-portrait. The differential equation governing the dynamics contains two
terms: a constant production rate given by $r_m$, and a degradation rate
$\gamma_m m$, which depends on the current mRNA count. The main plot shows each
of the components in the $m$ vs. $dm/dt$ plot. Since $r_m$ does not depend on
the current number of mRNA, it gives a straight production rate as a function of
$m$. The total degradation rate depends linearly on the mRNA copy number, giving
a line with slope $\gamma_m$. When the two components are equal (bot lines
crossing), we obtain the steady-state mRNA value $m_{ss}$. The bottom line shows
a qualitative schematic of the flow of the system towards this steady state. The
further $m$ is from $m_{ss}$, the faster it moves towards this point as
schematized by the arrows' size. (C) Example of mRNA dynamics for different
initial conditions. Over time all curves converge to the steady-state mRNA value
$m_{ss}=r_m/\gamma_m$. For this plot $\gamma_m = 1$ and $r_m/\gamma_m = 10$. The
[Python code
(`ch1_fig01C.py`)](https://github.com/mrazomej/phd/blob/master/src/chapter_01/code/ch1_fig01C.py)
used to generate part (C) of this figure can be found on the thesis [GitHub
repository](https://github.com/mrazomej/phd).](ch1_fig01){#fig:ch1_fig01
short-caption="Minimal model of gene expression"}

So far, our model assumes a simple constant transcription rate $r_m$; let us
expand this term a little further to include regulation by a transcriptional
repressor further down the road. We know that for a transcriptional event to
occur, the RNA polymerase (RNAP) must bind to the promoter region and undergo a
series of irreversible steps, such as opening the double helix to initiate the
DNA sequence's copying into mRNA [@Browning2004]. But before these irreversible
steps take place, there is a chance that the RNAP falls off the promoter. If we
assume these irreversible steps take place on a much longer timescale compared
to the initial binding and unbinding of the RNAP on the promoter, we can
separate the time scale and investigate them independently. In particular, we
can write that mRNA production happens at a rate
$$
\text{mRNA production} = r_m \cdot p_{\text{bound}},
\label{eq:mRNA_prod}
$$
where we split the original production term into two steps: $p_{\text{bound}}$,
the probability of finding an RNAP bound to the promoter, and $r_m$ which
captures all of the irreversible downstream steps that take place once the RNAP
is engaged in a transcriptional event. A way to think about it--relevant to what
I am doing right now as I type my thesis--is to think that the speed at which I
type this document has to do with two things: The probability of me being
actively working on these notes times the rate at which I type these notes once
I engage in the activity. The reason this separation makes sense is that we can
include the effect of the regulation by a transcriptional repressor as a
reduction of the time (the probability) that the RNAP can be bound to the
promoter. Furthermore, since we are assuming that the binding and unbinding of
the RNAP happen at a timescale much faster than the downstream events, we can
assume this binding reaction is in quasi-equilibrium, for which we can use the
powerful theoretical framework of statistical mechanics. Let us now delve into
the basics of this physical theory.

### The unreasonable effectiveness of unrealistic simplifications

In the preface of the textbook *Molecular Driving Forces* Dill and Bromberg
introduce the idea of Statistical Mechanics as *the unreasonable effectiveness
of unrealistic simplifications* [@Dill2010]. Although one could make the case
that all of physics follows this description, it is undoubtedly evident that
statistical mechanics is a vivid example of how simple ideas can have profound
consequences. Statistical mechanics can be defined as the theory that, upon
assuming the atomic nature of matter, explains the phenomenology that classical
thermodynamics established from the interactions of the microscopic components
of a system [@Dill2010]. As with any other physical theory, statistical
mechanics is built from a set of *empirical* facts that define "axioms" that we
take to be true. In other words, as Feynman famously described to us: if we want
to come up with a new law of nature, there is a simple recipe that we must
follow:

1. We guess the law. Literally. The most profound understanding of our physical
   reality we have comes from educated guesses made after a careful observation
   of nature.

2. We compute the consequences of such a guess. That is why mathematical
   theories allow us to sharpen our statements about how we think nature works.

3. We compare with experiments/observations. The scientific revolution came
   about when, after the dark ages, we finally learned it was okay to say "we
   don't know."

In such a simple statement, Feynman tells us, lies the key to science
[@Feynman1965]. For our purpose of understanding the basis of statistical
mechanics, we will argue that Boltzmann's law gives the main law upon which the
field is founded
$$
\frac{P(E_1)}{P(E_2)} = \frac{e^{-E_1 / k_BT}}{e^{-E_2 / k_BT}}.
\label{eq:boltzmann_law}
$$
Let us unpack this equation. The main idea behind statistical mechanics is that
macroscopic observables (temperature and pressure in classic examples) are
emergent properties dictated by the dynamics of the system's microscopic
components. What Boltzmann's law tells us is that the relative probability of a
system in thermal equilibrium to be found in a particular microstate with energy
$E_1$ compared to being in a microstate with energy $E_2$ is given by an
exponential function of the negative energy of such microstate relative to the
thermal energy $k_BT$. The minus sign in the exponent comes from the fact that
states with negative energies are more favorable by convention in physics. Thus,
having a large negative energy has a high probability when taking the
exponential of minus such negative number. To provide concrete examples of what
a microstate can look like, [@Fig:ch1_fig02](A) shows three molecular systems
relevant to biology. In the first example, we have the classic ligand-receptor
binding problem; here, we imagine a solution can be discretized in space into a
series of small boxes. In each of these boxes, one and only one ligand molecule
can fit in. In principle, we can list all possible spatial arrangements of
ligands. We could then calculate the relative likelihood of finding the system
in any configurations as long as we can assign an energy value to each of them.
The second example focuses on ligand-gated ion channels. In this particular
system, we care about the ion channel's state--either open or closed--and the
ligands' binding configuration. If the channel responds to the ligand's
concentration by changing its probability of gating, we can calculate using
equilibrium statistical mechanics. Finally, the third example shows different
configurations of a small patch of the cell membrane. All deformations of a
membrane have energetic costs associated with them. By listing all possible
membrane configurations, we can calculate the most likely shape of a membrane
given the forces and stresses acting on it. 

The macroscopic states that we observe can then be thought of as a
coarse-graining of many microstates into a single macrostate. For example, in
the ligand-receptor binding case, we rarely would care about the specific
position of all the ligand molecules in the solution. What we would be
interested in is whether or not the ligand is bound to the receptor. We can
therefore define as our "macrostate" the particular configuration of the
receptor as schematically shown in [@Fig:ch1_fig02](B).

![**Boltzmann's law and the definition of a micro and macrostate.** (A) Top
panel: ligand-receptor binding microstates. Middle panel: ligand-gated ion
channel microstates. Bottom panel: membrane patch deformations. (B) Schematic of
the definition of a "macrostate." In the ligand-receptor binding problem, we
ignore all ligand molecules' spatial configuration and focus on the receptor's
binding state.](ch1_fig02){#fig:ch1_fig02 short-caption="Boltzmann's law and the
definition of a micro and macrostate"}

If we want to know the likelihood of finding a particular system in any specific
configuration, Boltzmann's law (Eq. $\ref{eq:boltzmann_law}$) is then telling us
a protocol we must follow: 

1. Enumerate all possible microstates in which the system can be found.

2. Compute the energy of each of these microstates.

3. Define the "macrostate" we care about by grouping all microstates that belong
   to the same energy.

4. Compute the Boltzmann factor. This factor, sometimes called the Boltzmann
   weight, is defined as the exponential of the negative energy divided by the
   thermal energy, as indicated in Eq. $\ref{eq:boltzmann_law}$.

To see this protocol in action, let us apply it to the calculation of
$p_{\text{bound}}$, the probability of finding an RNAP bound to the promoter. We
will go through each of the protocol steps and build up the "unrealistic
simplifications" that will allow us to make this calculation.

**1. Enumerate possible microstates.** We begin by making a drastic
coarse-graining of the bacterial genome. For us, a genome is simply made out of
boxes where the RNAP can bind. We imagine that there is a single site where RNAP
can bind specifically--the promoter of interest. There are also $N_{NS} \approx
5\times 10^6$ non-specific binding sites, one per basepair (bp) in the genome.
This means that because of the sequence-dependent interactions between the RNAP
molecule, and the DNA, the energy associated with specific binding to the gene
promoter is more favorable than the rest of the genome. We ignore the fact that
the RNAP footprint where it binds to the genome is roughly 30 bp. This
assumption is valid if the number of available RNAP molecules is much smaller
than the number of non-specific binding sites since it is improbable that two
RNAPs would fall next to each other by pure chance. A useful analogy for this
point is to think about sitting $\sim \text{few}\times 10$ people on a large
stadium with $\sim 10^4$ seats. If the seats are chosen randomly, we do not need
to worry about doing the sampling "without replacement" because the chances of
two people ending up with the same seat number are negligible. We also ignore
the possibility of RNAP not being bound to the genome. This assumption is
supported by experimental evidence on a particular type of *E. coli* mutant that
sheds lipid vesicles without segregating DNA into such vesicles. Mass
spectrometry analysis on these "min-cells" has shown that there are no RNAP
molecules to be found, implying that RNAPs are bound to DNA most if not all of
the time [@Bintu2005]. The exercise then consists of randomly choosing one box
for each of the $P$ polymerases available to bind. [@Fig:ch1_fig03] shows in the
first column two possible configurations of our coarse-grained genome.

**2. Compute the energy for each microstate.** Let us analyze the case where all
$P$ RNAP molecules are bound non-specifically to the genome. For simplicity, we
assume that RNAP binds to all $N_{NS}$ non-specific binding sites with the same
affinity. We assign this energy to be $\varepsilon_P^{(NS)}$. This assumption
could be relaxed and we could assign instead a distribution of non-specific
binding energies, as explored in [@Phillips2019]. But for now, we don't have to
worry about this complication. For the statistical mechanics' protocol the
assignment of binding energies does not come from some quantum first-principled
calculation or anything similar. We label the interaction of the RNAP and the
rest of the genome with a single value, $\varepsilon_P^{(NS)}$, that
coarse-grains all of the hydrogen bonds and other effects that go into this
physical process and gives an average energy. The calculation continues with
this "labeled energy," and, as we will see at the end, a very clean functional
form emerges. Since we have $P$ such polymerases bound non specifically, the
energy of any state with a similar configuration is then $P
\varepsilon_P^{(NS)}$ as shown in [@Fig:ch1_fig03] second column, top row.

**3. Define the "macrostate" we care about.** In a sense, when we speak about
macrostate, it does not necessarily mean something that we can macroscopically
observe. What it means is that we group a bunch of states that we take to be
functionally equivalent, as shown in [@Fig:ch1_fig02](B). In our case, we only
care about whether or not the RNAP is bound to our promoter of interest. The
configuration of the rest of the background sites is irrelevant to our question.
What this means in practice is that we must compute the degeneracy or
multiplicity of our state. In other words, for the *specific* state shown in the
first column/top row of [@Fig:ch1_fig03] we know its Boltzmann weight. Eq.
$\ref{eq:boltzmann_law}$ tells us that the probability of this particular
configuration takes the form
$$
P_{\text{state}} \propto e^{-\beta P \varepsilon_P^{(NS)}},
$$
where we define $\beta \equiv (k_BT)^{-1}$. The probability of this binding
configuration takes this form since the $P$ RNAP molecules are bound non
specifically. But every single arrangement in which all RNAPs are bound
non-specifically has the same Boltzmann weight. The question then becomes: in
how many of such microstates can the system exist? This is a combinatorics
question of the form: in how many different ways can I arrange $P$ molecules
into $N_{NS}$ boxes? Which of course, the answer is
$$
\text{\# states with all RNAPs bound non-specifically} = 
\frac{N_{NS}!}{P!(N_{NS} - P)!},
$$
as shown in the third column of [@Fig:ch1_fig03]. This multiplicity can be
simplified if we consider that $N_{NS} \gg P$. To more easily visualize how to
simplify this let us for a second assume $N_{NS} = 100$ and $P = 3$. Given the
definition of factorials this means that
$$
\frac{N_{NS}!}{(N_{NS} - P)!} = 
\frac{100\cdot 99\cdot 98\cdots97\cdots 2\cdot 1}{97\cdots2\cdot 1} = 
100\cdot 99\cdot 98.
$$
Given this result, we can simply state that $100\cdot 99\cdot 98 \approx 100^3$,
only making a three percent error ($100\cdot 99\cdot 98 / 100^3$ \approx 0.97).
Imagine $N_{NS}$ is in the order of $10^6$, then the error would become
negligible. That is why, as shown in the third column of [@Fig:ch1_fig03], we
can approximate
$$
\frac{N_{NS}!}{P!(N_{NS} - P)!} \approx \frac{N_{NS}^P}{P!}, \;
\text{for }N_{NS} \gg P.
$$
For our other "macrostate" we have the case where only one out of the $P$ RNAPs
is bound specifically for the promoter. We define the energy of this single RNAP
specifically binding to the promoter as $\varepsilon_{P}^{(S)}$. We assume that
the other $P-1$ RNAPs are bound non-specifically with the usual energy
$\varepsilon_{P}^{(NS)}$. The way to realize this state is then given by
$$
\small
\text{\# states with one RNAP bound specifically} = 
\frac{N_{NS}!}{(P - 1)!(N_{NS} - (P - 1))!} \approx
\frac{N_{NS}^{P-1}}{(P-1)!}.
$$
What these Boltzmann weights mean is that for us *any* state on which a single
RNAP is bound to the promoter while the rest are bound non specifically is
equivalent. Therefore the probability of finding the promoter occupied by an
RNAP would be of the form
$$
p_{\text{bound}} \propto e^{-\beta \epsilon_1} + e^{-\beta \epsilon_2} +
e^{-\beta \epsilon_3} + \cdots
$$
where $\epsilon_i$ is the energy of the $i^{\text{th}}$ state that has a single
RNAP bound to the promoter. But we established that all of the $\epsilon_i$
energies are the same. So instead of writing this long sum, we multiply the 
Boltzmann weight of a single state by the number of states with equivalent
energy, i.e., we multiply it by the state's multiplicity or degeneracy. The 
same logic applies for the states where none of the RNAPs are specifically bound
to the promoter.

**4. Compute the Boltzmann Factor.** The last step in the protocol is to follow
the recipe indicated by Eq. $\ref{eq:boltzmann_law}$. We exponentiate the
energy, with the caveat we mentioned on the last point that this time we
multiply by the multiplicity that we just computed. This is because we are
lumping together all microstates into a single functional macrostate. So the
Boltzmann weight for the unbound $\rho_{\text{unbound}}$ macrostate is given by
$$
\rho_{\text{unbound}} = \frac{N_{NS}^P}{P!}
e^{-\beta P \varepsilon_P^{(NS)}}.
$$
For the bound state, we have
$$
\rho_{\text{bound}} = \frac{N_{NS}^{P-1}}{(P-1)!}
e^{-\beta \left(\varepsilon_P^{(S)} +  (P - 1) \varepsilon_P^{(NS)}\right)}.
$$
For reasons that will become clear later in this chapter once we work with the
entropy and derive the Boltzmann distribution, we know that to compute the
probability of a specific microstate (or a macrostate), we simply take the
Boltzmann weight of the microstate and divide by the *sum* of all of the other
Boltzmann weights of the states available to the system. This sum of Boltzmann
weights place a very special role in statistical mechanics, and it is known as
the *partition function* of the system. Therefore, to calculate
$p_{\text{bound}}$ we compute
$$
p_{\text{bound}} = 
\frac{\rho_{\text{bound}}}{\rho_{\text{unbound}} + \rho_{\text{bound}}}.
$$
Substituting the Boltzmann weights we derived, we find
$$
p_{\text{bound}} = 
\frac{
\frac{N_{NS}^{P-1}}{(P-1)!}
e^{-\beta \left(\varepsilon_P^{(S)} +  (P - 1) \varepsilon_P^{(NS)}\right)}
}{
\frac{N_{NS}^{P-1}}{(P-1)!}
e^{-\beta \left(\varepsilon_P^{(S)} +  (P - 1) \varepsilon_P^{(NS)}\right)}
+ 
\frac{N_{NS}^P}{P!}
e^{-\beta P \varepsilon_P^{(NS)}}
},
$$
an algebraic nightmare. We can simplify this expression enormously by
multiplying the numerator and denominator by $\rho_{\text{unbound}}^{-1}$. Upon
simplification, we find the neat expression
$$
p_{\text{bound}} = 
\frac{
    \frac{P}{N_{NS}} e^{-\beta \Delta \varepsilon_P}
}{
    1 + \frac{P}{N_{NS}} e^{-\beta \Delta \varepsilon_P}
},
\label{eq:pbound_unreg}
$$
where $\Delta\varepsilon_P \equiv \varepsilon_P^{(S)} - \varepsilon_P^{(NS)}$.
This simple expression, known as the Langmuir isothermal binding curve, tells us
that the more RNAPs available (larger $P$), or the stronger the promoter is
(more negative $\Delta\varepsilon_P$), the more likely it is to find the
promoter bound by an RNAP, and according to Eq. $\ref{eq:mRNA_prod}$, the higher
the mRNA production. In the next section, we connect this model to experimental
measurements. 

![**Statistical Mechanics protocol for RNAP binding.** On a discretized genome
we follow the statistical mechanics' protocol to compute the Boltzmann weight of
each of the relevant microstates. The $P$ available RNAPs are assumed to have
two binding configurations: One specific binding to the promoter of interest
(with energy $\varepsilon_P^{(S)}$) and non-specific to any of the $N_{NS}$
non-specific binding sites (with energy
$\varepsilon_P^{(NS)}$).](ch1_fig03){#fig:ch1_fig03 short-caption="Statistical
Mechanics protocol for RNAP binding"}

### Figure 1 theory in gene regulation

We began this section with a simple model for the dynamics of mRNA production
and degradation. We then expanded our model to deconvolve the production term
into the rate at which mRNA is produced by RNAP, and the probability of finding
such RNAP bound to the promoter. To calculate this probability, we used the
statistical mechanics' protocol, which culminated in Eq.
$\ref{eq:pbound_unreg}$. So far, we are missing two important steps in our
logical construction that will lead us to specific quantitative predictions that
we can test experimentally:

1. The inclusion of a regulatory scheme via a transcriptional repressor.

2. The connection of the model with experimentally accessible quantities.

As hinted at earlier, for a transcriptional repressor, we imagine that the
repressor's effect on the regulation of the gene acts only through changes in
$p_{\text{bound}}$. To include the regulation, we add a series of microstates.
Rather than having only $P$ RNAP molecules to bind the genome, we also have $R$
repressors that can bind specifically and non-specifically. Through the same
statistical mechanics' protocol as for the previous case, we can arrive at the
Boltzmann weights shown for the three "macrostates" in [@Fig:ch1_fig04](A). For
the regulated case, we have that the probability of the promoter being bound by
an RNAP takes the form
$$
p_{\text{bound}}(R > 0) = 
\frac{
    \frac{P}{N_{NS}} e^{-\beta \Delta \varepsilon_P}
}{
    1 
    + \frac{P}{N_{NS}} e^{-\beta \Delta \varepsilon_P}
    + \frac{R}{N_{NS}} e^{-\beta \Delta \varepsilon_R}
},
\label{eq:pbound_reg}
$$
where $\Delta\varepsilon_R$ is the binding energy difference between the
repressor binding to a specific binding site and a non-specific one. Although
exciting and insightful, the quantities we have derived so far do not have an
immediate **quantitative** prediction we can connect with experimental
measurements. For example, for the regulated case, the steady-state mRNA count
takes the form
$$
m_{ss}(R > 0) = 
\frac{r_m}{\gamma_m} p_{\text{bound}}(R > 0) = 
\frac{r_m}{\gamma_m}
\frac{
    \frac{P}{N_{NS}} e^{-\beta \Delta \varepsilon_P}
}{
    1 
    + \frac{P}{N_{NS}} e^{-\beta \Delta \varepsilon_P}
    + \frac{R}{N_{NS}} e^{-\beta \Delta \varepsilon_R}
}.
$$
Determining $r_m$ or $\gamma_m$ directly from experiments, although possible,
represents an enormous technical challenge. A convenient metric we can use
instead is what we call the fold-change in gene expression. [@Fig:ch1_fig04](B)
shows a schematic representation of what we mean by the fold-change. This
ratiometric quantity normalizes the expression level of a gene with regulation
given by a transcriptional repressor by the expression level of the same gene in
the absence of the regulation--via a knock-out of the repressor gene, for
example. Mathematically this is defined as
$$
\text{fold-change} \equiv \frac{m_{ss}(R > 0)}{m_{ss}(R = 0)}.
$$
This expression is convenient because upon taking the ratio of these
steady-state mRNA counts, the ratio $r_m / \gamma_m$ drops out of the equation.
All we are left with is then the ratio of the $p_{\text{bound}}$s
$$
\text{fold-change} = \frac{p_{\text{bound}}(R > 0)}{p_{\text{bound}}(R = 0)}.
$$
Substitutin Eqs. $\ref{eq:pbound_unreg}$ and $\ref{eq:pbound_reg}$ results in
$$
\text{fold-change} = 
\frac{
    1 
    + \frac{P}{N_{NS}} e^{-\beta \Delta \varepsilon_P}
}{
    1 
    + \frac{P}{N_{NS}} e^{-\beta \Delta \varepsilon_P}
    + \frac{R}{N_{NS}} e^{-\beta \Delta \varepsilon_R}
}.
$$
We appeal to some experimental understanding of the bacterial proteome
composition [@Bremer1996;@Schmidt2016;@Grigorova2006]. RNAP copy number in *E.
coli* is of the order $P \sim 10^3-10^4$ [@Grigorova2006]. The binding affinity
of these promoters is of the order $\Delta\varepsilon_P \sim -2\pm 1\;k_BT$
[@Bintu2005]. Along with the value of $N_{NS}\sim 10^6$ This results in
$$
\frac{P}{N_{NS}} e^{-\beta \Delta \varepsilon_P} \approx
\frac{10^3}{10^6}e^{2.3} \approx
\frac{10^3 \cdot 10}{10^6} \approx 10^{-2} \ll 1,
$$
the so-called weak-promoter approximation For the repressor we have that most
repressors in *E. coli* are in the order of $R \sim 10$ [@Schmidt2016]. Their
binding affinities take values between $\Delta\varepsilon_R \sim -15 \pm 5\;
k_BT$ [@Bintu2005]. These numerical values then give
$$
\frac{R}{N_{NS}} e^{-\beta \Delta \varepsilon_R} \approx
\frac{10}{10^6}e^{15} \approx
\frac{10 \cdot 10^6}{10^6} \approx 10.
$$
If we implement these approximations, we can justify simplifying the fold-change
equation to take the form
$$
\text{fold-change} \approx
\left(
    1 + \frac{R}{N_{NS}} e^{-\beta\Delta\varepsilon_R}
\right)^{-1}.
\label{eq:fc}
$$
As shown in [@Fig:ch1_fig04](C), this expression points directly at two
experimental knobs that we can tune using molecular biology. We can modify the
number of repressors by changing the ribosomal binding site sequence (RBS) of
the repressor gene [@Garcia2011c]. What that means is that with a
sequence-dependent manner, the ribosome translates mRNAs according to a specific
region of the gene known as the RBS [@Chen2013]. Furthermore, we can change the
repressor's affinity for its binding site by mutating the binding site itself
[@Garcia2011c]. [@Fig:ch1_fig04](D) shows predictions of Eq. $\ref{eq:fc}$ for
different binding energies.

The model and the predictions presented here were worked out by Garcia and
Phillips in a classic publication in 2011 [@Garcia2011c]. In the next chapter we
build upon this theoretical scaffold to expand the predictive power of the model
by including the allosteric nature of the transcription factor that allows the
cells to change their genetic program upon the presence of an external molecule
as a response to the environment.

![**Figure 1 theory in transcriptional regulation.** (A) States and (normalized)
weights for the simple repression motif. The promoter can be found in three
states: 1) empty, 2) bound by an RNAP, 3) bound by a repressor. The same
statistical mechanics' protocol as in [@Fig:ch1_fig03] can be used to derive the
weights. (B) Schematic of the experimental determination of the fold-change in
gene expression. The expression level of a regulated strain is normalized by the
expression level of a strain with a repressor's knock-out. (C)
Experimentally accessible knobs predicted from the theoretical model. The number
of transcription factors can be tuned by changing the amount of protein produced
per mRNA. The binding energy of the repressor can be tuned by mutating the
basepairs in the binding site. (D) Fold-change as a function of the repressor
copy number for different binding energies. The [Python code
(`ch1_fig04D.py`)](https://github.com/mrazomej/phd/blob/master/src/chapter_01/code/ch1_fig04D.py)
used to generate part (C) of this figure can be found on the thesis [GitHub
repository](https://github.com/mrazomej/phd)](ch1_fig04){#fig:ch1_fig04
short-caption="Figure 1 theory in transcriptional regulation"}

### The second secret of life

TBD. Allostery

### All cells are equal, but some are more equal than others

One of the great discoveries that came from the single-cell biology revolution
where we began to measure individual cellular behavior rather than bulk
observations, was the discovery of the intrinsic cell-to-cell variability in
many aspects of biology, gene expression being the canonical example
[@Eldar2010]. This means that two cells with the same genome exposed to the same
conditions will not express the same number of mRNAs and proteins of any
specific gene. From a statistical physics perspective, this is not entirely
"surprising" since we know that a system can be found in many different
microstates as described in [@Fig:ch1_fig02](A). What is different here is that
a cell does not have an Avogadro number of mRNA (or, for that matter of
anything) in it, making these fluctuations more relevant. If we think of
fluctuations scaling as $\sqrt{N}$, that means that for an $N$ of $\approx$ ten
molecules or so, these variations can be significant in terms of the downstream
cellular behavior. Cells have to cope with these physical limitations on
precision, many times generating systems to actively buffer as much of the
"noise" as possible [@Voliotis2014a], other times using this intrinsic
variability to their advantage [@Balaban2004].

The central assumption behind the thermodynamic models of gene regulation that
we studied in the last section is that the gene expression is proportional to
the probability of finding an RNAP bound to the promoter
[@Gerland2002;@Bintu2005]. A consequence of this construction is that the
probability space--the set of all possible events captured by the
distribution--only looks at the state of the promoter itself, not at the state
of the mRNA copy number. That is why thermodynamic models of this kind do not
speak to the intrinsic cell-to-cell variability. For this, we need to use the
so-called chemical master equation framework [@Sanchez2013]. There are two ways
of thinking about the chemical master equation:

1. The "particle" point of view.

2. The occupation number point of view.

Depending on the context, we might want to use either of these approaches to
write down the master equation for our problem of interest. Let us look into
these two different ways of interpreting a master equation using our example of
a cell producing mRNA. For the particle point of view, schematized in
[@Fig:ch1_fig05](A), we imagine following the mRNA copy number $m$ of a single
cell. The number of mRNAs in the cell change stochastically from time point to
time point. On the one hand, there can be a transcriptional event that increases
the number of mRNAs, and on the other hand, an mRNA can be degraded, decreasing
the number of mRNAs. Suppose we imagine tracking this cell for a very long time.
In that case, we can quantify the fraction of the time that the cell spent with
zero mRNAs, one, two, and so on and from that, build the probability
distribution $P(m, t)$ of having $m$ mRNA at time $t$ (there is a subtle point
here of the process being memoryless, but I do not want to get into it). The
occupation number point of view, schematized in [@Fig:ch1_fig05](B), takes a
different perspective. For this, we imagine tracking not one but many cells
simultaneously. Each cell can either produce or degrade an mRNA on a short time
window, changing its total individual count. The probability $P(m, t)$ is then
built from counting how many cells out of the total have $m$ mRNAs.

Regardless of how we think about the chemical master equation, both of these
perspectives describe a Markov process. These are stochastic processes in which
a system transitions between different states, but the transitions between such
states are only governed by the transition rates between the states and the
current state of the system. In other words, a Markov process keeps no track of
the states it previously visited; the only factor that determines where is the
system going to head is its current state, and the transition rates out of such
state--that is why these are considered memoryless processes.
[@Fig:ch1_fig05](C) shows a schematic of what a Markov process looks like. The
schematic of the unregulated promoter indicates that there are two possible
reactions: an mRNA production with rate $r_m$ and degradation with rate
$\gamma_m$. The Markov process for this simple model can then be represented as
a series of nodes (representing the mRNA counts) connected with bi-directional
arrows (representing the transition rates between states) indicating that the
transitions can only take place between contiguous states.

In practice, the way we write down a chemical master equation is by a process
christened by Professor Jane Kondev as "spread-the-butter." The idea of spread
the butter is that some probability mass (the analogous of the butter) is to be
spread over the range of possible values (the equivalent of the toast) where
probability mass migrates in and out of a particular bin keeping the total
amount of probability to add up to one. The best way to explain this concept is
by following the schematic in [@Fig:ch1_fig05](D) and going through the math.
Let us imagine we are keeping track of a particular mRNA value $m$--the chemical
master equations are in reality, a system of many coupled equations, one for
each mRNA count. We want to write down an equation that describes what is the
probability of finding a cell with this particular count a small time window
into the future $P(m, t + \Delta t)$, where $t$ represents the time "right now,"
and $\Delta t$ is a tiny time increment. The master equation is nothing more
than a checks and balances notebook to keep track of all the flow of probability
mass in and out of the bin we are interested in, as shown in
[@Fig:ch1_fig05](D). Informally we would write the equation as
$$
P(m, t + \Delta t) = P(m, t)
+ \sum \left({\text{transitions from} \atop m'\text{ to }m}\right)
- \sum \left({\text{transitions from} \atop m\text{ to }m'}\right),
\label{eq:master_intuition}
$$
where we are describing the three main components that go into the equation for
$P(m, t + \Delta t)$:

1. The probability of having $m$ mRNA right now,

2. the inflow of probability from other copy numbers $m'$ via production and
   degradation,

3. the outflow of probability from $m$ to other copy numbers via production and
   degradation.

Taking our time window $\Delta t$ to be sufficiently small, we can focus only on
the two contiguous mRNA counts $m-1$ and $m+1$, and ignore the rest since jumps
from further counts become increasingly improbable as the time step gets
smaller. [@Fig:ch1_fig05](D) shows the four in- and outflows that can happen.
Let us rewrite Eq. $\ref{eq:master_intuition}$ following this schematic. If a
cell has $m - 1$ mRNA and during the time window $\Delta t$ produces one
molecule, then it passes from state $m - 1$ to state $m$. This transition
contributes to the inflow of probability mass by a factor $(r_m \Delta t) P(m-1,
t)$, where we can think of $r_m \Delta t$ as the probability of the
transcription event taking place during the time window, and this multiplies the
probability of having $m - 1$ mRNA to begin with. A similar argument can be made
for all transitions in and out of $m$ depicted in [@Fig:ch1_fig05](D), with the
only difference that as in Eq. $\ref{eq:m_t_Delta_t}$, the degradation of an
mRNA molecule is proportional to the total number of molecules. The resulting
equation for $P(m, t + \Delta t)$ then takes the form
$$
\begin{split}
P(m, t + \Delta t) = &P(m, t)
+ \overbrace{(r_m \Delta t) P(m-1, t)}^{m-1 \rightarrow m}
+ \overbrace{(\gamma_m \Delta t) (m + 1) P(m + 1, t)}^{m+1 \rightarrow m}\\
&- \overbrace{(r_m \Delta t) P(m, t)}^{m \rightarrow m+1}
- \overbrace{(\gamma_m \Delta t) m P(m, t)}^{m \rightarrow m-1}
\end{split}.
$$
We send the first term on the right-hand side to the left, divide both sides by
$\Delta t$ and take the limit when $\Delta t \rightarrow 0$. This gives us the
master equation we were searching for
$$
\frac{dP(m, t)}{dt} = 
r_m P(m - 1, t) + \gamma_m (m + 1) P(m + 1, t)
- r_m P(m, t) + \gamma_m m P(m, t).
\label{eq:master_simple}
$$

![**Chemical master equation in gene regulation.** (A-B) Different points of
view to understand the chemical master equation. (A) From the "particle" point
of view, we imagine following the time trajectory of *a single cell*. The
probability $P(m, t)$ of finding a cell with $m$ mRNAs at time $t$ is then
proportional to the time this cell spent with this number of molecules. (B) On
the occupation number point of view we imagine observing a large number of
isogenic cells (different colors represent the individuality of each cell). The
probability $P(m,t)$ is then interpreted as the fraction of the cells
representing such copy number exactly at time $t$. (C) Chemical master equations
mathematize the idea of Markov processes. For the case of the unregulated
promoter, the Markov process consists of a connection of an infinite number of
discrete states that cells can transition between by producing or degrading
mRNAs. (D) Spread-the-butter idea. Since probability is conserved, the central
bar's height changes slightly by having in- and outflow of probability mass from
the contiguous bins. The [Python code
(`ch1_fig05A.py`)](https://github.com/mrazomej/phd/blob/master/src/chapter_01/code/ch1_fig05A.py)
used to generate the plot in part (A) of this figure can be found on the thesis
[GitHub repository](https://github.com/mrazomej/phd).](ch1_fig05){#fig:ch1_fig05
short-caption="Chemical master equation in gene regulation"}
