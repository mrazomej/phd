## Entropy, information, and the math behind the bit

Central to the endeavor undertaken in this thesis is the idea that cells can
process information from the environment to up or down-regulate their genes to
generate an appropriate response to these external signals. Information as a
concept is a very plastic term that we commonly use to explain having helpful
knowledge to use to our advantage. Phrases such as "*that person carries so much
information in her brain. She truly knows everything!*" point at this somewhat
imprecise concept of what we mean by information. 

In 1948, while working at Bell Labs, Claude Shannon shocked the world with his
seminal work that would go to define the field of information theory
[@Shannon1948]. In his paper, Shannon gave us a precise mathematical definition
of information. To understand Shannon's logic better, we need to put it in the
context that he was thinking about: communication systems such as the telephone
or the telegraph. Although seemingly unrelated to our problem of cells sensing
the environment, these systems are incredibly powerful in their conceptual and
explanatory reach. For Shannon, the main problem of communication consisted of
reproducing a message emitted at one point in space and time with fidelity at a
different point. Usually, these messages carry with them *meaning* (otherwise,
why would we even want to send such messages) by which we typically mean that
the message "refers to or is correlated according to some system with certain
physical or conceptual entities" [@Shannon1948]. But for the task of engineering
a reliable communication system, this meaning is irrelevant--in the same way
that whatever the cell decides to do with the meaning of the signals obtained
from the environment can be thought as irrelevant for the biophysics of how the
signal is sensed.

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

[@Fig:ch1_fig07](B) shows an analogous schematic to [@Fig:ch1_fig07](A) with the
relevant components involved in the gene expression context that we focus on
this thesis. In our bacterial gene regulation model, the information source role
is played by a small molecule's environmental concentration. It is this signal
that the cells are trying to measure and respond to by up-regulating the
expression of a gene. This signal transmitter is the allosteric transcription
factor whose conformation depends on the concentration of the small molecule.
The receiver of the signal is the DNA promoter that orchestrates the protein
expression, which plays the receiver's role.

![**Abstract communication system.** (A) Reproduced from Shannon's original
seminal work [@Shannon1948]. The schematic shows an abstract communication
system with all the components. (B) Adaptation of the Shannon communication
system to the context of bacterial gene expression regulated by an allosteric
transcription factor.](ch1_fig07){#fig:ch1_fig07 short-caption="Abstract
communication system"}

Having this setup in mind, the question becomes: how do we mathematically define
what information is? This brings a somewhat subtle difference between two
related terms that many time are incorrectly used interchangeably: *Entropy* and
*Information*. Information allows the entity that possesses it to make
predictions with accuracy better than random, while entropy is a quantification
of how much we do not know [@Adami2016]. From these definitions, we see that
having information, therefore, reduces our uncertainty, i.e., reduces the
entropy. This means that for Shannon, the amount of information we have from a
source is related to that source's statistical structure and how much we can
predict the message produced by the source given our knowledge of this
statistical structure. Let us look at a concrete example: English text. We know
that written and spoken language is not completely random. For a message to be
meaningful, the choice of words has to come from a statistical structure that
obeys the language's grammar rules. The choice of letters within a word also
follows a certain statistical structure. Imagine the text shown in
[@Fig:ch1_fig08](A). This is arguably one of the most important and most
beautiful and insightful pieces of prose ever put together by a human mind as it
is the last paragraph of *On the Origin of Species* by Darwin. If we ignore the
paragraph's message and just quantify how often we find each of the 26 letters
in the english alphabet, we obtain a distribution like the one shown in
[@Fig:ch1_fig08](B). In this paragraph, we can see that the most common vowel is
*e* exactly as in English writ-large. This distribution $P(x)$ is therefore not
maximally random. In other words, if we were to put all letters in the paragraph
in a hat and pick one letter at random, we could bet more money on the outcome
being a letter *e* and make money over time given this knowledge of the
structure of the distribution. A maximally random distribution would be if all
letters appeared equally frequent in the paragraph, such that betting on any
letter coming out of the hat would give us equal chances of guessing right. If
instead of looking at the distribution of individual letters, we look at pairs
of letters, the distribution $P(x, y)$ over the paragraph is shown in
[@Fig:ch1_fig08](C). Here we can see that just as the letters were not
completely random, the letters' pairs are also not random. For example, if we
take the first letter of the pair to be *t*, we see that it is more commonly
followed by letter *h*. What this implies is that knowing that the first letter
of the pair was *t* reduced our uncertainty of what character could come next.
We would then say that knowing the first letter gave us *information* about the
possible outcomes of the second letter. In the next section we will follow
Shannon's original derivation to mathematically define both entropy and
information.

![**The statistical structure of the English language.** (A) Last paragraph of
*On the Origin of Species* by Charles Darwin. This serves as a rather nice
not-random text example. (B) Marginal distribution $P(x)$ of all 23 letters and
the space. The size of the squares is proportional to how often each letter
appears in the paragraph. (C) Joint distribution of pairs of characters $P(x,
y)$. All pairs of characters in (A) were counted to build this histogram. The
x-axis shows the first letter while the y-axis shows the second. For simplicity
in (B) and (C) all punctuation was ignored.](ch1_fig08){#fig:ch1_fig08
short-caption="The statistical structure of the English language"}

### Choice, Uncertainty, and Entropy

So far our discussion about what entropy and information mean has been vague and
not rigorous. In order to derive a formula to quantify these concepts we need to
get more mathematical. Let us assume that an information source (See
[@Fig:ch1_fig07](A)) produces elements of a message following a distribution
$\mathbf{p} = \{p_1, p_2, \ldots, p_n \}$, where each $p_i$ is the probability
of the $i^{\text{th}}$ element. These elements could be letters, words,
sentences, concentrations of a small molecule, etc. of which we have $n$
possibilities. What we are looking for is a metric $H(\mathbf{p})$ that
quantifies how much "choice" is involved in the selection of each of the
elements of the message. In other words, how uncertain we are about the message
that the information source will produce? We demand of our desired quantity
$H(\mathbf{p})$ that it satisfies three reasonable conditions [@Shannon1948]:

1. $H$ should be continuous in the $p_i$s. Different information sources might
   have slightly different distributions $\mathbf{p}$, nevertheless $H$ should
   still apply to all possible information sources.

2. If all of the elements of the distribution are equally likely, i.e. $p_i =
   1/n$, then $H$ should be a monotonic increasing function of $n$. What this
   means is that the more options to chose from, the more uncertain we are about
   the possible outcome. For example, we are more uncertain about the outcome of
   a fair 6-sided die than of a fair coin just because of the number of possible
   outcomes from each of these "information sources."

3. If the act of choosing one of the possible $n$ elements of our information
   source can be break down into two successive choices, the original $H$ should
   be the weighted sum of the individual $H$s. What thi means is illustrated in
   [@Fig:ch1_fig09](A)where we imagine having an information source with only
   $n=3$ choices, each with probabilities $\mathbf{p} = \{ 1/2, 1/3, 1/6\}$,
   which gives $H(1/2, 1/3, 1/6)$ for the left case. For the right case we
   imagine first choosing between the upper and the lower path, and then, if the
   lower path is chosen, a second choice is made. This property then demands
   that
   $$
    \overbrace{H(1/2, 1/3, 1/6)}^{\text{single choice}} = 
    \overbrace{H(1/2, 1/2)}^{\text{first choice}} +
    \overbrace{\frac{1}{2} H(1/3, 1/6)}^{\text{second choice}}.
   $$
   Another way to think about this property is that we want our metric of
   uncertainty $H$ to be *additive*.

We will now proof that the only functional form that satisfies all these three
properties is given by
$$
    H(\mathbf{p}) = - K \sum_{i=1}^n p_i \log p_i,
$$
where $K$ is a constant having to do with the units (choice of the logarithm
base). In order to proof this we will follow Shannon's original work. Just as
Shannon did, we imagine the problem of encoding a message. Let us for example
imagine encoding a message from the english alphabet into Morse code as
schematically depicted in [@Fig:ch1_fig09](B). In there we take letters in
english alphabet (*SOS*), run it through an encoding function $E$ and obtain the
message (...---...). This process of encoding can be though as taking a message
$m_x$ written in an alphabet $\mathcal{X} = \{x_1, x_2, \ldots, x_n \}$, (where
$n$ is 26 for the english alphabet) and converting it into a message $m_y$
written in a different alphabet $\mathcal{Y} = \{y_1, y_2, \ldots, y_m \}$
(where $m=2$ for Morse code since we only have dots and dashes). The encoding
function $E: \mathcal{X}^r \rightarrow \mathcal{Y}^t$ takes a message of length
$r$ (for our example $r=3$ with three letters, *SOS*) and translates it into a
message of size $t$ (in our example $t=9$) such that we then have
$$
m_y = E(m_x),
$$
i.e., the function $E$ takes messages in english alphabet as input and spits out
a message in Morse code. It is obvious then that the larger the message $m_x$ we
want to encode, the larger the corresponding message $m_y$ will be. Therefore we
have that
$$
L(m_y) \propto L(m_x),
$$
where $L(\cdot)$ is a function that counts the number of characters in a
message. We have then $L(m_x) = r$ and $L(m_y) = t$. Let us call $k$ this
proportionality constant such that
$$
L(m_y) = k L(m_x).
\label{eq:length_fn}
$$
The engineering problem of designing an efficient coding scheme is then reduced
to find the smallest $K$ possible.

The number of messages of size $r$ that can be encoded with the alphabet
$\mathcal{X}$ is given by $n^r$--because we have $m$ possible options to chose
from for each of the $r$ characters, resulting in $n\cdot n\cdot n\cdots = n^r$.
Likewise, the number of messages of size $t$ encoded with alphabet $\mathcal{Y}$
is $m^t$. We then demand for our coding scheme that the number of messages that
we can encode is at least the number of messages we could potentially send. In
other words, for our coding scheme to be able to take *any* message of size $r$
it must be true that the number of possible encoded messages is at least as 
large as the number of possible messages to encode. This is expressed as
$$
n^r \leq m^t.
$$
If our encoding did not satisfy this, we would have to increase $t$, i.e. the
number of characters that we use to encode our message. On the other extreme we
could come up with a ridiculously long encoding scheme (imagine a version of
Morse code where every letter is represented by 1000 dots and dashes). To avoid
this absurd scheme we bound the size of the encoded message to be as long as
necessary to encode all potential messages, but not any longer. This is
expressed as
$$
m^{t-1} < n^r \leq m^t.
\label{eq:ineq_messages}
$$
Let us now take the logarithm on our previous inequality--this preserves the
inequalities since $\log$ is a monotonically increasing function--finding
$$
(t - 1) \log(m) < r \log(n) \leq t \log(m).
$$
We are free to take the logarithm in any base it is convenient, therefore let
us use base $m$ for this, obtaining
$$
t - 1 < r \log_m(n) \leq t.
\label{eq:ineq_logm}
$$
Dividing Eq. $\ref{eq:ineq_logm}$ by $r$ gives
$$
\frac{t-1}{r} < \log_m(n) \leq \frac{t}{r}.
\label{eq:t_over_r}
$$
Let us stare at Eq. $\ref{eq:t_over_r}$. In Eq. $\ref{eq:ineq_messages}$ We
established $t$ as the minimum number of characters from alphabet $\mathcal{Y}$
needed to encode a message of length $r$ written with alphabet $\mathcal{X}$
characters (such as *SOS* turned into ...---...). This means That $t/r$ is the
average number of characters from alphabet $\mathcal{Y}$ per character from
alphabet $\mathcal{Y}$, i.e., the proportionality constant $k$ from Eq.
$\ref{eq:length_fn}$. This means that Eq. $\ref{eq:t_over_r$ implies
$$
\log_m(n) \leq k.
$$
In words, this means that a lower bound for the number of characters from
alphabet $\mathcal{Y}$ needed to encode a character from alphabet $\mathcal{X}$
is given by $\log_m(n)$. Furthermore, Eq. $\ref{eq:t_over_r}$ implies that
$$
\frac{t}{r} - \log_m(n) < \frac{t}{r} - \frac{(t-1)}{r},
\label{eq:t_over_r_diff}
$$
given that $(t-1)/r < \log_m(n)$. Simplifying Eq. $\ref{eq:t_over_r_diff}$
results in
$$
\frac{t}{r} - \log_m(n) < \frac{1}{r}
\Rightarrow k - \log_m(n) < \frac{1}{r}.
$$
We can therefore make $k$, the number of encoding characters

![**TBD.** .](ch1_fig09){#fig:ch1_fig09
short-caption="TBD"}
