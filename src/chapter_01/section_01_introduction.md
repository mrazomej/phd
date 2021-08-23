## Introduction

In his classic 1944 book *What is Life?*, Schrödinger brought to the attention
of the scientific community what he thought were two of the biggest challenges
we had ahead of us if we were to understand living systems in the same way we
understand the electromagnetic field or the universal law of gravitation
[@Schrodinger1992]. The idea that living organisms could be "accounted for" by
physics and chemistry brought with it a new agenda on what needed to be done to
transition from a qualitative and descriptive study of the phenomena of life to
a quantitative and predictive science in the spirit of the physical sciences.
Since the publication of the book, there has been an enormous amount of progress
on our understanding of living systems from a first-principles perspective,
nevertheless, 75 years later Schrödinger questions are still as relevant and as
vibrant as ever before [@Phillips2021]. 

One of the defining features of living organisms at all scales is their capacity
of gathering information from the environment, encode an internal representation
of the state of the environment, and generate a response based on this
information processing capacity. Researchers in the field of origins-of-life had
gone as far as declaring that life emerged when chemical systems underwent a
radical phase transition after which they were able to process and use
information and free energy [@Cronin2016]. So, although speculative, it is
highly probable that the physical theory fulfilling Schrödinger's vision of
accounting for the phenomena of life will be the physics of systems capable of
processing information [@Davies2019].

In this context, information does not take the generic concept of possessing
practical knowledge about something. In this thesis, we use a precise
mathematical definition of information [@Adami2016]. This formal definition
makes information a metric worth quantifying and predicting in various
biological contexts as theoretical studies suggest that natural selection might
act on the ability of an organism to process information [@Taylor2007]. Thus,
working out the physical details of how it is that organisms sense the
environment---this is, gather information about the state of the environment,
encode such information in some shape or form within their physical boundaries,
and take action based on this information---is at the core of the
state-of-the-art research in biophysics [@Bialek2012].

The present thesis is an effort towards this vision of understanding biological
systems as information processing machines. Our object of study will be gene
regulation in bacteria. This particular system has been the subject of study for
microbiologists and molecular biologists for decades, and we have come to learn
a lot about the microscopic mechanistic details of how bacteria turn on and off
their transcriptional machinery [@Browning2004]. In particular, we will focus on
what we think of as the "hydrogen atom" of gene regulation---the so-called
simple-repression motif (more on that in the next section). In physics, calling
something the hydrogen atom of $X$ means that for the area of study $X$, this
"something" represents a system simple enough to be amenable to analytical
models that standard mathematical methods can solve, but rich enough to capture
the general features of the phenomena. This simple genetic circuit will allow us
to write tractable mathematical models to guide our experimental efforts with
the ultimate goal of testing our understanding of such systems when predicting
how much information a bacterium can gather from the environment using this
genetic module.

Professional biophysicists might wish to skip the rest of this chapter as we
will lay the foundations needed for the rest of our enterprise. We will
introduce the basics of gene expression modeling and the mathematical concept of
information and work through every single physical and mathematical prerequisite
needed for the rest of the thesis. The following chapters are structured as
follows: Chapter 2 builds on a decade of understanding this hydrogen atom of
gene regulation and expands our models' predictive ability by including the
effect of environmental effectors. This means that we will consider how gene
regulation is affected by the presence of an extracellular inducer molecule.
Chapter 3 will expand even further our predictive capacities by building a model
capable of making predictions about the cell-to-cell variability inherent to all
signaling systems working at the molecular scale. Chapter 4 serves as a
Supporting Information section for Chapter 2, detailing every calculation and
every inference. Likewise, Chapter 5 expands on Chapter 3, explaining every
technical detail.
