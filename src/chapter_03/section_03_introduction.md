## Introduction

As living organisms thrive in a given environment, they are faced with constant
changes in their surroundings. From abiotic conditions such as temperature
fluctuations or changes in osmotic pressure, to biological interactions such as
cell-to-cell communication in a tissue or a bacterial biofilm, living organisms
of all types sense and respond to external signals. [@Fig:ch3_fig01](A) shows a
schematic of this process for a bacterial cell sensing a concentration of an
extracellular chemical. At the molecular level, where signal transduction
unfolds mechanistically, there are physical constraints on the accuracy and
precision of these responses given by intrinsic stochastic fluctuations
[@Nemenman2010]. This means that two genetically identical cells exposed to the
same stimulus will not have identical responses [@Eldar2010].

One implication of this noise in biological systems is that cells do not have an
infinite resolution to distinguish signals. Consequently, there is a one-to-many
mapping between inputs and outputs. Furthermore, given the limited number of
possible outputs, there are overlapping responses between different inputs. This
scenario can be mapped to a Bayesian inference problem where cells try to infer
the state of the environment from their phenotypic response, as schematized in
[@Fig:ch3_fig01](B). The question then becomes this: how can one analyze this
probabilistic, rather than deterministic, relationship between inputs and
outputs? The abstract answer to this question was worked out in 1948 by Claude
Shannon who, in his seminal work, founded the field of information theory
[@Shannon1948]. Shannon developed a general framework for how to analyze
information transmission through noisy communication channels. In his work,
Shannon showed that the only quantity that satisfies three reasonable axioms for
a measure of uncertainty was of the same functional form as the thermodynamic
entropy--thereby christening his metric the information entropy [@MacKay2003].
Based on this information entropy, he also defined the relationship between
inputs and outputs known as the mutual information. The mutual information $I$
between input $c$ and output $p$, given by
$$
I = \sum_c P(c) \sum_p P(p \mid c) \log_2 \frac{P(p \mid c)}{P(p)},
\label{eq:ch3_eq01}
$$
quantifies how much we learn about the state of the input $c$ given that we get
to observe the output $p$. In other words, the mutual information can be thought
of as a generalized correlation coefficient that quantifies the degree to which
the uncertainty about a random event decreases given the knowledge of the
average outcome of another random event [@Kinney2010].

It is natural to conceive of scenarios in which living organisms can better
resolve signals might have an evolutionary benefit, making it more likely that
their offspring will have a fitness advantage [@Taylor2007]. In recent years
there has been a growing interest in understanding the theoretical limits on
cellular information processing [@Bialek2005; @Gregor2007], and in quantifying
how close evolution has pushed cellular signaling pathways to these theoretical
limits [@Tkacik2008; @Dubuis2013; @Petkova2019]. While these studies have
treated the signaling pathway as a "black box," explicitly ignoring all the
molecular interactions taking place in them, other studies have explored the
role that molecular players and regulatory architectures have on these
information processing tasks [@Rieckh2014; @Ziv2007; @Voliotis2014a;
@Tostevin2009; @Tkacik2011; @Tkacik2008a; @Tabbaa2014]. Despite the great
advances in our understanding of the information processing capabilities of
molecular mechanisms, the field still lacks a rigorous experimental test of
these detailed models with precision measurements on a simple system in which
physical parameters can be perturbed. In this work, we approach this task with a
system that is both theoretically and experimentally tractable in which
molecular parameters can be varied in a controlled manner.

Over the last decade, the dialogue between theory and experiments in gene
regulation has led to the predictive power of models not only over the mean
level of gene expression but the noise as a function of relevant parameters such
as regulatory protein copy numbers, the affinity of these proteins to the DNA
promoter, as well as the extracellular concentrations of inducer molecules
[@Golding2005; @Garcia2011c; @Vilar2013; @Xu2015]. These models based on
equilibrium and non-equilibrium statistical physics have reached a predictive
accuracy level such that, for simple cases, it is now possible to design
input-output functions [@Brewster2012; @Barnes2019]. This opens the opportunity
to exploit these predictive models to tackle how much information genetic
circuits can process. This question lies at the heart of understanding the
precision of the cellular response to environmental signals. [@Fig:ch3_fig01](C)
schematizes a scenario in which two bacterial strains respond with different
levels of precision to three possible environmental states, i.e., inducer
concentrations. The overlap between the three different responses precisely
determines the resolution with which cells can distinguish different inputs.
This is analogous to how the point spread function limits the ability to resolve
two light-emitting point sources.

In this work, we follow the same philosophy of theory-experiment dialogue used
to determine model parameters to predict from first principles the effect that
biophysical parameters such as transcription factor copy number and protein-DNA
affinity have on the information processing capacity of a simple genetic
circuit. Specifically, to predict the mutual information between an
extracellular chemical signal (input $c$, isopropyl
$\beta$-D-1-thiogalactopyranoside or IPTG in our experimental system) and the
corresponding cellular response in the form of protein expression (output $p$),
we must compute the input-output function $P(p \mid c)$. To do so, we use a
master-equation-based model to construct the protein copy number distribution as
a function of an extracellular inducer concentration for different combinations
of transcription factor copy numbers and binding sites. Having these
input-output distributions allow us to compute the mutual information $I$
between inputs and outputs for any arbitrary input distribution $P(c)$. We opt
to compute the channel capacity, i.e., the maximum information that can be
processed by this gene regulatory architecture, defined as [Eq:ch3_eq01]
maximized over all possible input distributions $P(c)$. By doing so we examine
the physical limits of what cells can do in terms of information processing by
harboring these genetic circuits. Nevertheless, given the generality of the
input-output function $P(p \mid c)$ we derive, the model presented here can be
used to compute the mutual information for any arbitrary input distribution
$P(c)$. All parameters used for our model were inferred from a series of studies
that span several experimental techniques [@Garcia2011c; @Jones2014a;
@Brewster2014; @Razo-Mejia2018], allowing us to make parameter-free predictions
of this information processing capacity [@Phillips2019].

These predictions are then contrasted with experimental data, where the channel
capacity is inferred from single-cell fluorescence distributions taken at
different inducer concentrations for cells with previously characterized
biophysical parameters [@Garcia2011c; @Razo-Mejia2018]. We find that our
parameter-free predictions quantitatively track the experimental data up to a
systematic deviation. The lack of numerical agreement between our model and the
experimental data poses new challenges towards having a foundational,
first-principles understanding of the physics of cellular decision-making.

<!-- The reminder of the paper is organized as follows. In we define the minimal
theoretical model and parameter inference for a simple repression genetic
circuit. discusses how all parameters for the minimal model are determined from
published datasets that explore different aspects of the simple repression
motif. computes the moments of the mRNA and protein distributions from this
minimal model. In we explore the consequences of variability in gene copy number
during the cell cycle. In this section we compare experimental and theoretical
quantities related to the moments of the distribution, specifically the
predictions for the fold-change in gene expression (mean expression relative to
an unregulated promoter) and the gene expression noise (standard deviation over
mean). follows with reconstruction of the full mRNA and protein distribution
from the moments using the maximum entropy principle. Finally uses the
distributions from to compute the maximum amount of information that the genetic
circuit can process. Here we again contrast our zero-parameter fit predictions
with experimental inferences of the channel capacity. -->

![**Cellular signaling systems sense the environment with different degrees of
precision**. (A) Schematic representation of a cell as a noisy communication
channel. From an environmental input (inducer molecule concentration) to a
phenotypic output (protein expression level), cellular signaling systems can be
modeled as noisy communication channels. (B) We treat cellular response to an
external stimulus as a Bayesian inference of the state of the environment. As
the phenotype (protein level) serves as the internal representation of the
environmental state (inducer concentration), the probability of a cell being in
a specific environment given this internal representation $P(c \mid p)$ is a
function of the probability of the response given that environmental state $P(p
\mid c)$. (C) The precision of the inference of the environmental state depends
on how well cells can resolve different inputs. For three different input levels
(left panel), the green strain responds more precisely than the purple strain
since the output distributions overlap less (middle panel). This allows the
green strain to make a more precise inference of the environmental state given a
phenotypic response (right panel).](ch3_fig01){#fig:ch3_fig01
short-caption="Cellular signaling systems sense the environment with different
degrees of precision"}
