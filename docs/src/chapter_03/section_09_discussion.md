## Discussion 

Building on Shannon's formulation of information theory, there have been
significant efforts using this theoretical framework to understand the
information processing capabilities of biological systems, and the evolutionary
consequences for organisms harboring signal transduction systems
[@Bergstrom2004; @Taylor2007; @Tkacik2008; @Polani2009; @Nemenman2010;
@Rivoire2011]. Recently, with the mechanistic dissection of molecular signaling
pathways, significant progress has been made on the question of the physical
limits of cellular detection and the role that features such as feedback loops
play in this task [@Bialek2005; @Libby2007; @Tkacik2011; @Rhee2012a;
@Voliotis2014a]. But the field still lacks a rigorous experimental test of these
ideas with precision measurements on a system that is tractable both
experimentally and theoretically.

In this chapter, we take advantage of the recent progress on the quantitative
modeling of input-output functions of genetic circuits to build a minimal model
of the simple repression motif [@Phillips2019]. By combining a series of studies
on this circuit spanning diverse experimental methods for measuring gene
expression under a myriad of different conditions, for the first time, we
possess complete *a priori* parametric knowledge---allowing us to generate
parameter-free predictions for processes related to information processing. Some
of the model parameters for our kinetic formulation of the input-output function
are informed by inferences made from equilibrium models. We use the fact that if
both kinetic and thermodynamic languages describe the same system, the
predictions must be self-consistent. In other words, if the equilibrium model
can only make statements about the mean mRNA and mean protein copy number
because of the way these models are constructed, those predictions must be
equivalent to what the kinetic model has to say about these same quantities.
This condition, therefore, constrains the values that the kinetic rates in the
model can take. To test whether or not the equilibrium picture can reproduce the
predictions made by the kinetic model, we compare the experimental and
theoretical fold-change in protein copy number for a suite of biophysical
parameters and environmental conditions ([@Fig:ch3_fig03](C) upper row). The
agreement between theory and experiment demonstrates that these two frameworks
can indeed make consistent predictions.

The kinetic treatment of the system brings with it increasing predictive power
compared to the equilibrium picture. Under the kinetic formulation, the
predictions are not limited only to the mean but any of the moments of the mRNA
and protein distributions. Furthermore, our formulation in terms of dynamical
equations allows us to account for the time-varying nature of the moments of the
mRNA and protein copy numbers. Specifically, since the protein mean lifetime is
comparable with the cell cycle length, the protein copy number does not reach a
steady-state over the cell cycle duration. Accounting for this effect increases
the expected cell-to-cell variability when measuring non-synchronized cells. We
first test these novel predictions by comparing the noise in protein copy number
(standard deviation/mean) with experimental data. Our minimal model predicts the
noise up to a systematic deviation. The physical or biological origins of this
discrepancy remain an open question. In that way, the work presented here
exposes the status quo of our understanding of gene regulation in bacteria,
posing new questions to be answered with future model refinements. We then
extend our analysis to infer entire protein distributions at different input
signal concentrations using the maximum entropy principle. This means that we
compute moments of the protein distribution and then use these moments to build
an approximation to the full distribution. These predicted distributions are
then compared with experimental single-cell distributions, as shown in
[@Fig:ch3_fig04](B) and [Sec. 5.5](#sec:ch5_sec06). Again, although our minimal
model systematically underestimates the width of the distributions, it informs
how changes in parameters such as protein copy number or protein-DNA binding
affinity will affect the full probabilistic input-output function of the genetic
circuit to a multiplicative constant. We then use our model to predict the
information processing capacity.

By maximizing the mutual information between input signal concentration and
output protein distribution over all possible input distributions, we predict
the channel capacity of the system over a suite of biophysical parameters such
as varying repressor protein copy number and repressor-DNA binding affinity.
Although there is no reason to assume the simplified synthetic circuit we used
as an experimental model operates optimally given the distribution of inputs,
the relevance of the channel capacity comes from its interpretation as a metric
of the physical limit of how precise of an inference cells can make about what
the state of the environment is. Our model, despite the systematic deviations,
makes non-trivial predictions such as the existence of an optimal repressor copy
number for a given repressor-DNA binding energy, predicting the channel capacity
up to an additive constant (see [@Fig:ch3_fig05]). The origin of this optimal
combination of repressor copy number and binding energy differs from previous
publications in which an extra term associated with the cost of producing
protein was included in the model [@Tkacik2011]. This optimal parameter
combination is a direct consequence of the fact that the LacI repressor cannot
be fully deactivated [@Razo-Mejia2018]. This implies that as the number of
repressors increases, a significant number of them are still able to bind to the
promoter even at saturating concentrations of inducer. This causes all of the
input-output functions to be shifted towards low expression levels, regardless
of the inducer concentration, decreasing the amount of information that the
circuit can process. Interestingly, the number of bits predicted and measured in
our system is similar to that of the gap genes in the *Drosophila* embryo
[@Dubuis2013]. Although this is a suggestive numerical correspondence that sets
current experimental data on the information processing capacity of genetic
circuits between 1 and 2 bits, more work is required to fully understand the
effect that different regulatory architectures have on the ability to resolve
different signals.

We consider it important to highlight the limitations of the work presented
here. The previously discussed systematic deviation for the noise and skewness
of the predicted distributions (see [Sec. 5.8](#sec:ch5_sec09)), and therefore
of the predicted distributions and channel capacity, remains an unresolved
question. Our current best hypothesis for the origin of this unaccounted noise
pertains to cell-to-cell variability in the central dogma machinery. More
specifically, our model does not account for changes in RNAP and sigma factor
copy numbers, changes in ribosome numbers, and even the variability in the
repressor copy number. This possibility deserves to be addressed in further
iterations of our minimal model. Also, as first reported in [@Razo-Mejia2018],
our model fails to capture the steepness of the fold-change induction curve for
the weakest repressor binding site (see [@Fig:ch3_fig03](B)). Furthermore, the
minimal model in (A), despite being widely used, is an oversimplification of the
physical picture of how the transcriptional machinery works. The coarse-graining
of all the kinetic steps involved in transcription initiation into two effective
promoter states---active and inactive---ignores potential kinetic regulatory
mechanisms of intermediate states [@Scholes2017]. Moreover, it has been argued
that even though the mRNA count distribution does not follow a Poisson
distribution, this effect could be caused by unknown factors, not at the level
of transcriptional regulation [@Choubey2018].

The findings of this work open the opportunity to accurately test intriguing
ideas that connect Shannon's metric of how accurately a signaling system can
infer the state of the environment with Darwinian fitness [@Taylor2007].
Beautiful work along these lines has been done in the context of the
developmental program of the early *Drosophila* embryo [@Tkacik2008;
@Petkova2019]. These studies demonstrated that the input-output function of the
pair-rule genes works at channel capacity, suggesting that selection has acted
on these signaling pathways, pushing them to operate at the limit of what the
physics of these systems allow. Our system differs from the early embryo because
we have a tunable circuit with variable amounts of information processing
capabilities. Furthermore, compared with the fly embryo in which the organism
tunes both the input and output distributions over evolutionary time, we have
experimental control of the distribution of inputs that the cells are exposed
to. Consequently, this means that instead of seeing the final result of the
evolutionary process, we would be able to set different environmental challenges
and track over time the evolution of the population. These experiments could
shed light on the suggestive hypothesis of information bits as a trait on which
natural selection acts. We see this exciting direction as part of the overall
effort in quantitative biology of predicting evolution [@Lassig2017].
