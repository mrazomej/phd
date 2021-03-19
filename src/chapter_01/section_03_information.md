## Entropy, information, and the math behind the bit

Central to the endeavor undertaken in this thesis is the idea that cells can
process information from the environment in order to up- or down-regulate their
genes to generate an appropriate response to these external signals. Information
as a concept is a very plastic term that we commonly use to explain the
existence of having useful knowledge that we can then use in our advantage.
Phrases such as "*that person carries so much information in her brain. She
truly knows everything!*" point at this rather imprecise concept of what we mean
by information.

In 1948, while working at Bell labs, Claude Shannon shocked the world with his
seminal work that would go to define the field of information theory
[@Shannon1948]. In his paper Shannon gave us a precise mathematical definition
of information. To understand Shannon's logic better we need to put it in the
context that he was thinking about: communication systems such as the telephone
or the telegraph. These systems, although seemingly unrelated to our problem of
cells sensing the environment, are incredibly powerful in their conceptual and
explanatory reach. For Shannon the main problem of communication consisted in
reproducing a message emitted at one point in space and time with fidelity at
a different point. Usually these messages carry with them *meaning* (otherwise 
why would we even want to send such messages) by which we usually mean that the
message "refers to or is correlated according to to some system with certain 
physical or conceptual entities" [@Shannon1948]. But for the task of 
engineering a reliable communication system this meaning is irrelevant--in the
same way that whatever the cell decides to do with the meaning of the signals
obtained from the environment can be thought as irrelevant for the biophysics
of how the signal is sensed.

As shown schematically in [@Fig:ch1_fig07](A) from Shannon's original work, a
communication system essentially consists of five components:

1. An **information source** which produces a message (or sequence of messages)
   to be communicated to the receiving terminal.

2. A **transmitter** which takes the message, converts it into a suitable signal
   compatible with the communication channel.

3. The **channel** that is the medium used to transmit the signal from the
   transmitter to the receiver.

4. The **receiver** in charge of inverting the operation done by the 
   transmitter, reconstructing the original message.

5. The **destination** for whom the message is intended.

[@Fig:ch1_fig07](B) shows an analogous schematic to [@Fig:ch1_fig07](A) with
the relevant components involved in the gene expression context that we focus
on this thesis. In our bacterial gene regulation model the information source
role is played by the environmental concentration of a small molecule. It is 
this signal that the cells are trying to measure and respond to by up-regulating
the expression of a gene. The transmitter of this signal is the allosteric
transcription factor whose conformation depends on the concentration of the
small molecule. The receiver of the signal is the DNA promoter that orchestrates
the expression of the protein, which plays the role of the receiver.

![**Abstract communication system.** (A) Reproduced from Shannon's original
seminal work [@Shannon1948]. The schematic shows an abstract communication
system with all the components. (B) Adaptation of the Shannon communication
system to the context of bacterial gene expression regulated by an allosteric
transcription factor.](ch1_fig07){#fig:ch1_fig07 short-caption="Abstract
communication system"}

Having this setup in mind, the question becomes: how do we define mathematically
what information is? This brings a rather subtle difference between two related
terms that many time are incorrectly use interchangeably even in textbooks:
*Entropy* and *Information*. Information is that which allows the entity that
possess it to make predictions with accuracy better than random, while entropy
is a quantification of how much we do not know [@Adami2016]. From these
definitions we see that having information therefore reduces our uncertainty,
i.e., reduces the entropy.

For Shannon, the amount of information that
a source carries with it is related to the statistical structure of that source.
The more structure (less random) a source is, the more predictable it becomes,
and therefore the less we learn when we get to observe the outcome--consequently 
having less information. At first this can be counterintuitive; how can 
something structured and predictable carry little information? Books are really
structured, nevertheless we think of them as carrying a lot of information. But
what Shannon asks us to do is to think of information as something related to
how surprised we are of seeing certain outcome. In the example of the book, 
imagine the most boring book that consists of 400 pages of just one word 
repeated over and over again. Nobody can deny that book is incredibly structured
and completely predictable, but do we gain any information when we turn through
the pages reading the same word? Definitely not.


![**The randomness of the English language.** (A) Last paragraph of *On the
Origin of Species* by Charles Darwin. This serves as a rather nice not-random
text example. (B) Marginal distribution $P(x)$ of all 23 letters and the space.
The size of the squares is proportional to how often each letter appears in the
paragraph. (C) Joint distribution of pairs of characters $P(x, y)$. All pairs of
characters in (A) were counted to build this histogram. The x-axis shows the
first letter while the y-axis shows the second. For simplicity in (B) and (C)
all punctuation was ignored.](ch1_fig08){#fig:ch1_fig08 short-caption="The
randomness of the English language"}
