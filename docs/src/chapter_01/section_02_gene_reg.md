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

### The unreasonable effectiveness of unrealistic simplifications

On the preface of the textbook *Molecular Driving Forces* Dill and Bromberg
introduce the idea of Statistical Mechanics as *the unreasonable effectiveness
of unrealistic simplifications [@Dill2010]. Although one could make the case
that all of physics follows this description, it is certainly evident that
statistical mechanics is a vivid example of how simple ideas can have profound
consequences. Statistical mechanics can be simply defined as the theory that,
upon assuming the atomic nature of matter, explains the phenomenology that
classical thermodynamics established [@Dill2010] from the interactions of the
microscopic components of a system. As with any other physical theory,
statistical mechanics is built from a set of *empirical* facts that define
axioms that we take to be truth. In other words, as Feynman famously described
to us: if we want to come up with a new law of nature there is a simple recipe
that we must follow:
1. We guess the law. Literally. The deepest understanding of physical reality we
   have up to now come from educated guesses made after a careful observation of
   nature.
2. We compute the consequences of such guess. That is why mathematical theories
   allow us to sharpen our statements about how nature works
3. We compare with experiments/observations. The scientific revolution came
   about when, after the dark ages, we finally learned it was okay to say "we
   don't know."

In such simple statement, Feynman tells us, lies the key to science
[@Feynman1965]. For our purpose of understanding the basis of statistical
mechanics one could argue that the main law upon which the field is founded is
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
microstate can look like, [@Fig:ch1_fig01](A) shows three molecular systems
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
them. So listing all possible microstates of the system we can calculate what is
the most likely configuration that we could find a membrane given the forces and
stresses acting on the membrane.

The macroscopic states that we get to observe can then be thought of as a
coarse-graining of many microstates into a single macrostate. For example, in
the case of the ligand-receptor binding, we rarely would care about the specific
position of all the ligand molecules in the solution. What we would be
interested in is whether or not the ligand is bound to the receptor. We can
therefore define as our "macrostate" the particular configuration of the
receptor as schematically shown in [@Fig:ch1_fig01](B).

![**Boltzmann's law and the definition of a micro and macrostate.** (A) Top
panel: ligand-receptor binding microstates. Middle panel: ligand-gated ion
channel microstates. Bottom panel: membrane patch deformations. (B) Schematic of
the definition of a "macrostate." In the ligand-receptor binding problem we
ignore the spatial configuration of all ligand molecules, and focus on the
binding state of the receptor.](ch1_fig01){#fig:ch1_fig01
short-caption="Boltzmann's law and the definition of a micro and macrostate"}

If we want to know the likelihood of finding a particular system in any specific
configurationBoltzmann's law (Eq. $\ref{eq:boltzmann_law}$) is then telling us a
protocol: 
1. Enumerate all possible microstates in which the system can be found.
2. Compute the energy of each of these microstates.
3. Compute the Boltzmann factor by exponentiating minus the energy divided by
   the thermal energy.
