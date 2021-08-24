## Entropy, Information, and the Math Behind the Bit

Central to the endeavor undertaken in this thesis is the idea that cells can
process information from the environment to up or down-regulate their genes to
generate an appropriate response to these external signals. Information as a
concept is a very plastic term that we commonly use to explain having helpful
knowledge to use to our advantage. Phrases such as "*That person carries so much
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
a reliable communication system, this meaning is irrelevant---in the same way
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
relevant components involved in the gene expression context that we focus on in
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
predict the source's message given our knowledge of this statistical structure.
Let us look at a concrete example: English text. We know that written and spoken
language is not completely random. For a message to be meaningful, the choice of
words has to come from a statistical structure that obeys the language's grammar
rules. The choice of letters within a word also follows a certain statistical
structure. Let us look at the text shown in [@Fig:ch1_fig08](A). This is
arguably one of the most important and most beautiful pieces of prose ever put
together by a human mind as it is the last paragraph of *On the Origin of
Species* by Darwin. If we ignore the paragraph's message and just quantify how
often we find each of the 26 letters in the English alphabet, we obtain a
distribution like the one shown in [@Fig:ch1_fig08](B). This paragraph shows
that the most common vowel is *e*, exactly as in English writ-large. This
distribution $P(x)$ is therefore not maximally random. In other words, if we
were to put all letters in the paragraph in a hat and pick one letter at random,
we could bet more money on the outcome being a letter *e* and make money over
time given this knowledge of the structure of the distribution. 

A maximally random distribution would be if all letters appeared equally
frequent in the paragraph, such that betting on any letter coming out of the hat
would give us equal chances of guessing right. If instead of looking at the
distribution of individual letters, we look at pairs of letters, the
distribution $P(x, y)$ over the paragraph is shown in [@Fig:ch1_fig08](C). Here
we can see that, just as the letters were not completely random, the pairs of
letters are also not random. For example, if we take the first letter of the
pair to be *t*, we see that it is more commonly followed by the letter *h*. This
implies that knowing that the first letter of the pair was *t* reduced our
uncertainty of what character could come next. We would then say that knowing
the first letter gave us *information* about the possible outcomes of the second
letter. In the next section, we will follow Shannon's original derivation to
define both entropy and information mathematically.

![**The statistical structure of the English language.** (A) Last paragraph of
*On the Origin of Species* by Charles Darwin. This serves as a rather nice
not-random text example. (B) Marginal distribution $P(x)$ of all 26 letters and
space. The size of the squares is proportional to how often each letter appears
in the paragraph. (C) Joint distribution of pairs of characters $P(x, y)$. All
pairs of characters in (A) were counted to build this histogram. The x-axis
shows the first letter while the y-axis shows the second. For simplicity in (B)
and (C) all punctuation was ignored. The [Python code
(`ch1_fig08.py`)](https://github.com/mrazomej/phd/blob/master/src/chapter_01/code/ch1_fig08.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/mrazomej/phd).](ch1_fig08){#fig:ch1_fig08
short-caption="The statistical structure of the English language"}

### Choice, Uncertainty, and Entropy

So far, our discussion about what entropy and information mean has been vague
and not rigorous. To derive a formula to quantify these concepts, we need to get
more mathematical. Let us assume that an information source (see
[@Fig:ch1_fig07](A)) produces elements of a message following a distribution
$\mathbf{p} = \{p_1, p_2, \ldots, p_n \}$, where each $p_i$ is the probability
of the $i^{\text{th}}$ element. These elements could be letters, words,
sentences, basepairs, concentrations of a small molecule, etc., of which we have
$n$ possibilities. What we are looking for is a metric $H(\mathbf{p})$ that
quantifies how much "choice" is involved in the selection of each element of the
message. In other words, how uncertain we are about the message that the
information source will produce at random? We demand our desired quantity
$H(\mathbf{p})$ to satisfy three reasonable conditions [@Shannon1948]:

1. $H$ should be continuous in the $p_i$s. Different information sources might
   have slightly different distributions $\mathbf{p}$, nevertheless $H$ should
   still apply to all possible information sources.

2. If all of the elements of the distribution are equally likely, i.e., $p_i =
   1/n$, then $H$ should be a monotonic increasing function of $n$. This means
   that the more options to choose from, the more uncertain we are about the
   possible outcome. For example, we are more uncertain about the outcome of a
   fair 6-sided die than of a fair coin just because of the number of possible
   outcomes from each of these "information sources."

3. If the act of choosing one of the possible $n$ elements of our information
   source can be broken down into two successive choices, the original $H$
   should be the weighted sum of the individual $H$s. What this means is
   illustrated in [@Fig:ch1_fig09](A) where we imagine having an information
   source with $n=3$ choices, each with probabilities $\mathbf{p} = \{ 1/2, 1/3,
   1/6\}$, which gives $H(1/2, 1/3, 1/6)$ for the left case. For the right case,
   we imagine first choosing between the upper and the lower path, and then, if
   the lower path is chosen, a second choice is made. This property then demands
   that
   $$
    \overbrace{H(1/2, 1/3, 1/6)}^{\text{single choice}} = 
    \overbrace{H(1/2, 1/2)}^{\text{first choice}} +
    \overbrace{\frac{1}{2} H(1/3, 1/6)}^{\text{second choice}}.
   $$
   Another way to think about this property is that we want our metric of
   uncertainty $H$ to be *additive*.

![**Shannon's theorem.** (A) One of the properties of a reasonable metric for
uncertainty is that we can partition choices into multiple steps, and the
resulting uncertainty should remain the same. (B) Example of coding functions
$E$. The English alphabet can be converted into Morse code. Amino acids can be
encoded in codons. (C) Partitioning of $2^3$ equally likely choices into three
decision steps, each with two choices. Eight different amino acids can be
selected using two schemes: 1) each of the eight codons is chosen at random with
equally likely chances, or 2) the codon is built by choosing one basepair at the
time. (D) Partitioning of unequal choices. Given the redundancy of the genetic
code, for equally likely codons, the resulting amino acid has different
probabilities being chosen.](ch1_fig09){#fig:ch1_fig09 short-caption="Shannon's
theorem"}

We will now prove that the only functional form that satisfies all these three
properties is given by
$$
    H(\mathbf{p}) = - K \sum_{i=1}^n p_i \log p_i,
$$
where $K$ is a constant having to do with the units (choice of the logarithm
base). To prove this, we will follow Shannon's original work. We imagine the
problem of encoding a message. For example, imagine encoding a message from the
English alphabet into Morse code, or a protein sequence into the corresponding
mRNA sequence, as schematically depicted in [@Fig:ch1_fig09](B). In there, we
take letters in the English alphabet (*SOS* for the English alphabet, *MGF* for
the protein), run it through an encoding function $E$ and obtain the message
(...- - -... for the Morse code, *AUGGGCUUC* for the mRNA). This process of
encoding can be thought of as taking a message $m_x$ written in an alphabet
$\mathcal{X} = \{x_1, x_2, \ldots, x_n \}$, (where $n$ is 26 for the English
alphabet, and 20 for the number of amino acids) and converting it into a message
$m_y$ written in a different alphabet $\mathcal{Y} = \{y_1, y_2, \ldots, y_m \}$
(where $m=2$ for Morse code since we only have dots and dashes, and $m=4$ for
the mRNA with 4 possible nucleotides). The encoding function $E: \mathcal{X}^r
\rightarrow \mathcal{Y}^t$ takes a message of length $r$ (for our exmaples
$r=3$) and translates it into a message of size $t$ (in our examples $t=9$)
such that we then have
$$
m_y = E(m_x).
$$
Obviously, the larger the message $m_x$ we want to encode, the larger the
corresponding message $m_y$ will be. Therefore we have that
$$
L(m_y) \propto L(m_x),
\label{eq:length_proportionality}
$$
where $L(\cdot)$ is a function that counts the number of characters in a
message. An essential difference between both of the examples in
[@Fig:ch1_fig09](B) is that, for the English to Morse code case, the number of
dots and dashes for different letters is different (*e*$\rightarrow$.,
*x*$\rightarrow$-..-). Meanwhile, for the amino acid to codon case, every single
codon has the same length. Let us focus for now on this second coding scheme
where every character from alphabet $\mathcal{X}$ is encoded with the same
number of characters from alphabet $\mathcal{Y}$. We have then $L(m_x) = r$ and
$L(m_y) = t$. Let us call $k$ the proportionality constant from Eq.
$\ref{eq:length_proportionality}$ such that
$$
L(m_y) = k L(m_x).
\label{eq:length_fn}
$$
The number of messages of size $r$ that can be encoded with the alphabet
$\mathcal{X}$ is given by $n^r$---because we have $n$ possible options to chose
from for each of the $r$ characters, resulting in $n\cdot n\cdot n\cdots = n^r$.
Likewise, the number of messages of size $t$ encoded with alphabet $\mathcal{Y}$
is $m^t$. We then demand from our coding scheme that the number of messages we
can encode is at least the number of messages we could potentially send. In
other words, for our coding scheme to be able to take *any* message of size $r$,
it must be true that the number of possible encoded messages is at least as
large as the number of possible messages to encode. This demand is expressed as
$$
n^r \leq m^t.
$$
If our encoding did not satisfy this, we would have to increase $t$, i.e., the
number of characters we use to encode our message. For example, if codons were
made out of only two basepairs, the genetic code would not be able to code for
all 20 amino acids plus the stop codons. On the other extreme, we could develop
a ridiculously long encoding scheme (imagine a version of the genetic code where
1000 basepair represented a single amino acid). To avoid this absurd scheme, we
bound the encoded message's size to be as long as necessary to encode all
potential messages, but not any longer. This bound is expressed as
$$
m^{t-1} < n^r \leq m^t.
\label{eq:ineq_messages}
$$
Let us now take the logarithm on our previous inequality---this preserves the
inequalities since $\log$ is a monotonically increasing function---finding
$$
(t - 1) \log(m) < r \log(n) \leq t \log(m).
$$
We are free to choose the logarithm base we find convenient; therefore, let us
use base $m$ for this, obtaining
$$
t - 1 < r \log_m(n) \leq t.
\label{eq:ineq_logm}
$$
Dividing Eq. $\ref{eq:ineq_logm}$ by $r$ gives
$$
\frac{t-1}{r} < \log_m(n) \leq \frac{t}{r}.
\label{eq:t_over_r}
$$
Let us stare at Eq. $\ref{eq:t_over_r}$. In Eq. $\ref{eq:ineq_messages}$, we
established $t$ as the minimum number of characters from alphabet $\mathcal{Y}$
needed to encode a message of length $r$ written with alphabet $\mathcal{X}$
characters (such as *MGF* turned into *AUGGGCUUC* as in [@Fig:ch1_fig09](B)).
This means that, for the case where all symbols use the same number of
characters when encoded, $t/r$ is the number of characters from alphabet
$\mathcal{Y}$ per character from alphabet $\mathcal{X}$, i.e., the
proportionality constant $k$ from Eq. $\ref{eq:length_fn}$. This means that Eq.
$\ref{eq:t_over_r}$ implies
$$
\log_m(n) \leq k.
$$
In other words, a lower bound for the number of characters from alphabet
$\mathcal{Y}$ needed to encode a character from alphabet $\mathcal{X}$ is given
by $\log_m(n)$. For the amino acid to codon case, the minimum number of letters
in a codon would be $\log_4(20) \approx 2.16 > 2$. This shows why we could not
encode all 20 amino acids with two basepair long codons. Furthermore, Eq.
$\ref{eq:t_over_r}$ implies that
$$
\frac{t}{r} - \log_m(n) < \frac{t}{r} - \frac{(t-1)}{r},
\label{eq:t_over_r_diff}
$$
given that $(t-1)/r < \log_m(n)$. Simplifying Eq. $\ref{eq:t_over_r_diff}$
results in
$$
\frac{t}{r} - \log_m(n) < \frac{1}{r}
\Rightarrow k - \log_m(n) < \frac{1}{r}.
\label{eq:logn_1_over_r}
$$
Therefore, we can make $k$, the number of encoding characters, as arbitrarily
close to $\log_m(n)$ as we want by increasing the length of the message being
encoded, i.e., making $r \rightarrow \infty$. This would imply a genetic code,
not for individual amino acids but entire polypeptides. This scheme would not
work biologically; nevertheless, this mathematical limit will help us find the
functional form of our desired function $H(\mathbf{p})$. 

Coming back to the function $H$, let us define
$$
A(n) \equiv H \left(\frac{1}{n}, \frac{1}{n}, \frac{1}{n}, \ldots\right),
$$
as the maximum possible value of $H$ when all outcomes are equally likely.
Property 2 tells us that $A(n)$ increases monotonically with the length of the
message. This means that if we apply the function $A(\cdot)$ to the terms in Eq.
$\ref{eq:ineq_messages}$, we conserve the inequality, i.e.,
$$
A(m^{t-1}) < A(n^r) \leq A(m^t).
\label{eq:A_ineq}
$$
Using Property 3, we can divide the $n^r$ possible choices into $r$ independent
decisions, each with $n$ options to chose from. This property is depicted in
[@Fig:ch1_fig09](C). On the left, it shows we can choose from eight different
codons that code for $2^3 = 8$ different amino acids. On the right, we can
choose base by base, building up the codon in three consecutive decisions, each
with two equally likely choices, for a total of $2\cdot 2 \cdot 2 = 8$ possible
outcomes. This division of choices allows us to rewrite Eq. $\ref{eq:A_ineq}$ as
$$
(t - 1) A(m) < r A(n) \leq t A(m),
\label{eq:A_div}
$$
because of our requirement of the uncertainty $H$ being an additive property.
For the example in [@Fig:ch1_fig09](C), at each of the three decision steps, the
uncertainty is given by $A(2)$. Given that the uncertainty is additive, for each
of the routes, our total uncertainty is given by
$$
A(2) + A(2) + A(2) = 3 A(2),
$$
therefore $A(2^3) = 3 A(2)$. Dividing Eq. $\ref{eq:A_div}$ by $r$ results in
$$
\frac{(t - 1)}{r} A(m) < A(n) \leq \frac{t}{r} A(m).
$$
Since $\frac{(t - 1)}{r} A(m) < A(n)$, it is also true that
$$
\frac{t}{r} A(m) - A(n) < A(m) \left(\frac{t}{r} - \frac{(t-1)}{r} \right).
$$
Simplifying terms, we are left with
$$
\frac{t}{r} A(m) - A(n) < \frac{1}{r} A(m).
$$
Dividing both sides by $A(m)$, we find 
$$
k - \frac{A(n)}{A(m)} < \frac{1}{r}.
\label{eq:A_1_over_r}
$$
We can make the ratio $A(n)/A(m)$ as close to $k$ as we want by making $r$
larger. This equation looks shockingly similar to Eq. $\ref{eq:logn_1_over_r}$,
but what is the connection? On the one, hand Eq. $\ref{eq:logn_1_over_r}$ is the
result of imposing the condition that our coding scheme must be able to encode
any possible message from one alphabet $\mathcal{X}$ to another alphabet
$\mathcal{Y}$. This condition leads us to the conclusion that the number of
characters from alphabet $\mathcal{Y}$ needed to encode the characters from
alphabet $\mathcal{X}$ (the constant $k$) can be made as arbitrarily close to
$\log_m(n)$ as we want by writing a code, not for individual characters
(individual amino acids), but for sequences of characters (polypeptides). On the
other hand, Eq. $\ref{eq:A_1_over_r}$ is a direct consequence of the three
logical properties we imposed on our uncertainty metric $H$. These properties
led us to conclude that, whatever our uncertainty function for the equally
likely choices $A(\cdot)$ is, the ratio of the uncertainties for each of our two
alphabets $A(n)/A(m)$ approaches the same constant $k$ as we make the encoded
message longer. Since both $\log_m(n)$ and $A(n)/A(m)$ approach $k$ as $r$
grows, we can conclude that
$$
\frac{A(n)}{A(m)} \rightarrow \frac{\log_m(n)}{\log_m(m)}\; \text{ as }
r \rightarrow \infty.
$$
We wrote the ratio $\log_m(n)/\log_m(m)$ because our choice of the logarithm
base was arbitrary. Therefore, more generally, we have
$$
\frac{A(n)}{A(m)} \rightarrow \frac{\log(n)}{\log(m)}\; \text{ as }
r \rightarrow \infty,
$$
for any base. This convergence only takes place if and only if
$$
A(n) = K \log(n),
$$
where $K$ is some constant. This is quite beautiful. What we just demonstrated
is that the functional form for the uncertainty metric we are after scales as
the logarithm of the number of possible characters in our alphabet. We know that
our uncertainty function $H(1/n, 1/n, \ldots)$ is a function of $1/n$ rather
than of $n$. This is easily fixed by using the properties of logarithms, writing
$$
H\left(\frac{1}{n}, \frac{1}{n}, \ldots \right) = 
-K \log\left(\frac{1}{n} \right).
\label{eq:entropy_equally}
$$
The general form of Shannon's entropy is starting to show up. After all, for the
case where all choices are equally likely, we have $p_i = 1/n$. We can therefore
write
$$
H\left(\frac{1}{n}, \frac{1}{n}, \ldots \right) = -K 
\sum_{i=1}^n \frac{1}{n} \log\left(\frac{1}{n} \right).
$$

Let us generalize the proof for cases where choices are not equally likely. To
continue with the amino acid to codon encoding example, we now consider the
genetic code's redundancy. Given that there are $4^3 = 64$ possible codons,
multiple codons map to the same amino acid. An example of three amino acids that
share the first letter is depicted on [@Fig:ch1_fig09](D). The diagram on the
left shows a total of nine different codons; two of such codons code for
asparagine (*N*), three for isoleucine (*I*), and four for threonine (*T*). A
way to express the asymmetry between the choices is to have each codon as an
independent and equally likely choice, as depicted on the middle diagram of
[@Fig:ch1_fig09](D). Let us define the total number of codons 
$$
N = \sum_{i=1}^n n_i,
\label{eq:sum_amino}
$$
where $n_i$ counts the number of codons for amino acid $i$, and $n$ is the total
number of amino acid choices. Let us call $H_1$ the uncertainty of this set of
equal choices. From Eq. $\ref{eq:entropy_equally}$, we know that the resulting
uncertainty function $H_1$ is of the form 
$$
H_1 = K \log\left(\sum_{i=1}^n n_i \right) = K \log (N),
$$
since all codons are equally likely.

Although each codon is equally likely, the resulting amino acid is not. The
probability of amino acid *I* in this case is the number of codons encoding it
(two) divided by the total number of codons in the example (nine). In general,
we assume that each of the $n$ choices has a probability
$$
p_i = \frac{\text{\# codons for amino acid }i}{\text{total \# of codons}} =
\frac{n_i}{N}.
\label{eq:p_i_amino}
$$
By Property 3 of our function $H$, we can partition the codon's choice into two
consecutive decisions (not three since the first codon is the same for all amino
acids in this example). This partitioning is shown on the right diagram of
[@Fig:ch1_fig09](D). The uncertainty $H_2$ for this case has two contributions,
one for each of the decisions
$$
H_2 = 
\overbrace{H(p_1, p_2, \ldots, p_n)}^{\text{first choice}} + 
\overbrace{K \sum_{i=1}^n p_i \log n_i}^{\text{second choice}}.
$$
The first decision has an unknown functional form we are trying to figure out.
The second choice consists of choosing between $n_i$ equally likely bases for
the codon's last position, each weighted by the probability of going to this
particular branch (the one that defines the amino acid) as demanded by Property
3. But whether or not we choose each codon on a single decision or in two steps,
the uncertainty of this event is the same. This means that $H_1 = H_2$ as
Property 3 requires. This equality results in
$$
K \log (N) = H(p_1, p_2, \ldots, p_n) + K \sum_{i=1}^n p_i \log(n_i).
$$
Solving for $H(p_1, p_2, \ldots, p_n)$ results in
$$
H(p_1, p_2, \ldots, p_n) = K \left[ 
    \log N - \sum_{i=1}^n p_i \log(n_i)
\right].
$$
Using Eq. $\ref{eq:sum_amino}$ results in
$$
H(p_1, p_2, \ldots, p_n) = - K \left[ 
    \sum_{i=1}^n p_i \log(n_i)
    - \log\left( \sum_{i=1}^n n_i \right)
\right].
$$
Since probabilities must be normalized, i.e., $\sum_{i=1}^n p_i = 1$, we can
write
$$
H(p_1, p_2, \ldots, p_n) = - K \left[ 
    \sum_{i=1}^n p_i \log(n_i)
    - \sum_{i=1}^n p_i \log\left( \sum_{i=1}^n n_i \right)
\right].
$$
Using the property of logarithms, we can rewrite this as
$$
H(p_1, p_2, \ldots, p_n) = - K \left[ 
    \sum_{i=1}^n p_i 
    \log\left( \frac{n_i}{\sum_{i=1}^n n_i} \right)
\right].
$$
Using Eq. $\ref{eq:p_i_amino}$, we find the expected result
$$
H(p_1, p_2, \ldots, p_n) = - K \sum_{i=1}^n p_i \log p_i.
\label{eq:shannon_result}
$$

Let us dissect this result. We began this derivation by stating three logical
properties that a metric for uncertainty should have. The properties could be
summarized simply as 1) the function exists for all possible $p_i$s, 2) the
uncertainty grows as the number of possible outcomes grows, and 3) the
uncertainty must be additive. We thought about a coding scheme to encode a
message written in an alphabet into a different one. We demanded that our coding
scheme should be able to encode *any* message we want, and this led us to
conclude that the average number of characters needed to encode each character
on the original message can approach $\log_m(n)$, where $n$ is the number of
characters in the original alphabet and $m$ is the number of characters in the
encoding alphabet. We then used the properties of our desired uncertainty
function and found a non-obvious connection between the number of characters
needed to pass from one alphabet to another and the uncertainty on the message.
When we generalized this analysis to cases where not all outcomes are equally
likely, we arrived at Eq. $\ref{eq:shannon_result}$, the so-called Shannon
entropy. This is Shannon's theorem, and what it shows is that Eq.
$\ref{eq:shannon_result}$ is the only function that satisfies the three very
reasonable conditions we established for an uncertainty measurement.

To gain intuition on what this equation is telling us, let us look at two
examples. In our first example, we will think about the simplest random process:
a coin toss. To compute how unpredictable the outcome of our simple coin toss
is, we can use Eq. $\ref{eq:shannon_result}$. For this particular case, we only
have two possible outcomes---heads with probability $p$ or tails with
probability
$1 - p$. The resulting entropy is of the form
$$
H = - p \log(p) - (1 - p) \log(1 - p).
\label{eq:entropy_coin}
$$
[@Fig:ch1_fig10](A) plots Eq. $\ref{eq:entropy_coin}$ as a function of the
probability of heads $p$. Notice that the curve is concave with a minimum at
$p=0$ and $p=1$ and a maximum at $p=1/2$. This shape should make intuitive sense
given that Eq. $\ref{eq:entropy_coin}$ quantifies how unpredictable the outcome
of tossing the coin is. If the coin toss's outcome is always heads (p=1) or
always tails (p=0), there is no uncertainty about the resulting face. The more
both outcomes become (the closer $p$ gets to $1/2$), the more unpredictable the
random even is. One mathematical subtly here is that for $p=1$ or $p=0$, we have
to compute $0 \times \log(0)$, which is undefined. For these values of $p$ we
take $0 \times \log(0) = 0$ since the limit where $x \rightarrow 0^+$ converges
to zero. Notice that the units on the $y$-axis are given in bits. These units
mean that we used base two for our logarithms. An easy way to think about what a
bit means is the number of *yes*/*no* questions one would need to ask on average
to infer the random event's outcome. For a coin, all we need is a single
question (therefore one bit) to know what the outcome was.

For our second example, we go back to the mRNA steady-state distribution we
derived in Eq. $\ref{eq:mRNA_steady}$. We found that for our simple one-state
DNA promoter, the steady-state distribution resulted in Poisson with a mean
$\langle m \rangle = r_m / \gamma_m$. [@Fig:ch1_fig10](B) shows the entropy of
this Poisson distribution as a function of the mean mRNA. We see a quick initial
increase in this entropy up to $\langle m \rangle \approx 20$, after which there
is a much less steep increment. Imagine we sample a random cell from one of
these Poisson distributions. Using the interpretation of bits again as the
number of *yes/no* questions, what [@Fig:ch1_fig10](B) tells us is that if the
promoter produces $\approx 10$ mRNA on average, it will take on average 3.5 of
these questions to infer the number of mRNA for random cell. For an average of
$\approx 20$ mRNA, it would take four questions, and for an average of $\approx
60$ mRNA, five questions. These questions would be of the form "*is it greater
than the average?*" or "*is it less than or equal to 1/3 of the average?*," and
so on.

![**Shannon entropy in action.** (A) The entropy of a coin as a function of the
probability of heads $p$. The entropy is maximum when the coin is fair, i.e.,
$p=0.5$, meaning that this is the most unpredictable coin one could have. (B)
The entropy of the steady-state mRNA distribution as derived in Eq.
$\ref{eq:mRNA_steady}$ as a function of the mean mRNA copy number. The point
shows the entropy of the distribution shown in the inset. Bot figures use base 2
for the logarithm, resulting in units of bits for the entropy. The [Python code
(`ch1_fig10.py`)](https://github.com/mrazomej/phd/blob/master/src/chapter_01/code/ch1_fig10.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/mrazomej/phd).](ch1_fig10){#fig:ch1_fig10
short-caption="Shannon entropy in action"}

### Information Theory and Statistical Mechanics

Our result in Eq. $\ref{eq:shannon_result}$ is of the same functional form as
the thermodynamic entropy. The story goes that Shannon was discussing this
concept with his friend John von Neumann. It was von Neumann who allegedly
convinced Shannon of calling his metric of randomness *entropy* under the
argument that nobody understands the concept. But the fact that the functional
forms are the same is too suggestive to dismiss a potential connection between
these concepts immediately. It was until much later that E. T. Jaynes formalized
ways to link both ideas [@Jaynes1957]. Nevertheless, Jaynes himself strongly
discourages people from trying to map one concept to the other explicitly. In
his book "*Probability Theory: The Logic of Science*," Jaynes warns the reader
about failing to distinguish information entropy, which is a property of the
mathematical object we call a probability distribution, and the *experimental
entropy* of thermodynamics, which is instead a property of the state of the
system as defined by experimentally measurable quantities such as volume,
temperature, pressure, magnetization, etc. Jaynes goes on to say: "*they should
never have been called by the same name; the experimental entropy makes no
reference to any probability distribution, and the information entropy makes no
reference to thermodynamics*" [@Jaynes2003].

When Jaynes makes such strong remarks about the disconnection between both
entropy concepts, he strictly refers to the classical thermodynamic definition.
This classical definition of entropy, due to Clausius, refers to the inability
of any thermal engine to convert all of the input energy into useful work.
Clausius defined a new quantity $S$ as the amount of energy per unit temperature
unavailable to do work. To understand this idea is to realize that from the
energy liberated in gasoline combustion on a car engine, we only end up
extracting $\approx 20\%$ of the energy to move the car. The other $80\%$ is
lost into heating the engine and the environment. But this is not because the
engineers are using poor designs. The second law of thermodynamics on its
classical definition states that nothing in the universe can convert $100\%$ of
the energy into useful work; there will always be residual energy that gets
turned into heat.

At the time, the existence of atoms was not widely accepted by the scientific
community. But then came Boltzmann and the statistical mechanics' conceptual
revolution. The giant leap in our understanding of why the second law of
thermodynamics does not allow the total conversion of energy into useful work
came with Boltzmann's revolutionary entropy idea. Boltzmann hypothesized that
matter was made out of atoms. Therefore, everything we can observe and measure
macroscopically about any system results from the microscopic configuration of
all the atoms that make up the system. Furthermore, many microscopic
arrangements are indistinguishable at our macroscopic scale (recall the
microstate and macrostate concept in [@Fig:ch1_fig02]). This line of reasoning
led Boltzmann to the law we stated in Eq. $\ref{eq:boltzmann_law}$. This law and
all of the classic results from statistical mechanics are founded on several
assumptions about the microscopic scale processes' reversibility. In other
words, for Boltzmann's law to be "a legit law of nature," it must be the case
that if we play a movie featuring a single atom moving around the system, the
same movie played in reverse should be as equally likely to happen.

But it might be the case that the assumptions underlying statistical mechanics
laws are not the most fundamental constructs of reality. As we will show next,
we can derive a classic result of statistical mechanics from a completely
different premise having to do more with statistical inference rather than
physical laws of motion governing atoms. This becomes a circular argument where
some physicists have the laws of motion as the defining foundation on which to
base statistical mechanics laws is better. For others, having an
information-theoretic justification for statistical mechanics independent of the
underlying physical laws is more appealing. At the end of the day is a matter of
taste. Having said all of this, let us delve into the connection between
information-theoretic entropy and the Boltzmann distribution.

We already used the Boltzmann distribution when we computed the probability of
an RNAP molecule being bound to the promoter $p_{\text{bound}}$. The Boltzmann
distribution applies to systems in thermodynamic equilibrium in contact with a
heat bath at a constant temperature. Think of a small Eppendorf tube ($\approx
2$ mL) that we perfectly seal before submerging it into the ocean. The tube's
temperature will equilibrate with that of the ocean, but the ocean's temperature
will not be affected by the tube's presence. Submerging the tube into the
reservoir allows the total energy of the tube not to be fixed. Sometimes
the tube can borrow energy from the ocean; sometimes, it can give energy to it.
The Boltzmann distribution precisely dictates the likelihood of such energy
states. The probability of a state with energy $E_i$ is given by
$$
P(E_i) = \frac{e^{-\beta E_i}}{\mathcal{Z}},
$$
where, as before, $\beta \equiv (k_BT)^{-1}$. $\mathcal{Z}$ is the partition
function defined by the sum of the Boltzmann weight for all possible
microstates, i.e.,
$$
\mathcal{Z} \equiv \sum_{\text{states}} e^{-\beta E_i},
$$
where the sum is taken over all microstates available to the system. This
equation is equivalent to Eq. $\ref{eq:pbound_unreg}$ and Eq.
$\ref{eq:pbound_reg}$. We can derive this functional form from the so-called
maximum entropy principle. This framework is expanded more in Chapter 5 of this
thesis. But for our purposes here, the idea is that we are trying to make a
"best guess" of what a distribution looks like, given limited information. For
our Eppendorf tube inside the ocean, we are thinking about the distribution of
all of the molecules' microstates inside the tube. Experimentally, we never get
to observe any of the microstates of the system. But we know that the
probability of each microstate depends on its energy, as Boltzmann told us. Let
us say we can measure the average energy $\langle E \rangle$ of our little
Eppendorf tube. What is then the optimal guess of the functional form of the
distribution that does not use any information we do not have at hand? For
example, we cannot say that there is only one microstate available to the system
with energy $\langle E \rangle$, because that constrains the possibilities of
the system, and measuring the average energy does not lead to such a conclusion.
The next best case we can do is to maximize the Shannon entropy, subject to this
constraint on the average energy. This makes sense because, as we derived in the
previous section, the Shannon entropy is the only functional form that satisfies
our properties for a metric of uncertainty. Maximizing the Shannon entropy leads
then to a maximally uninformative distribution. Including the constraints when
implementing this maximization guarantees that we use all that we know about the
distribution and nothing else.

Given Property 1 of our function $H$, the Shannon entropy is continuous on the
individual probabilities' values $p_i$. This means that we can maximize the
Shannon entropy by taking its derivative with respect to $p_i$ and equating it
to zero. This operation does not include the constraints we have on the values
of the probabilities of each microstate. Let us say that each microstate
available to the system with energy $E_i$ has a probability $p_i$ of happening.
The constraint on the average energy is given by
$$
\langle E \rangle = \sum_{\text{states}} E_i p_i,
$$
where again, the sum is taken over all possible microstates. Furthermore, we
know that the probability distribution must be normalized. This means that
$$
\sum_{\text{states}} p_i = 1.
$$
To include these constraints in our optimization, we can use the Lagrange
multipliers technique. We refer the reader to any introductory text on
multivariate calculus for a quick refresher of this technique. We proceed by
defining a Lagrangian $\mathcal{L}$ of the form
$$
\mathcal{L}(p_1, p_2,\ldots, p_N, \beta, \mu) =
\overbrace{- \sum_{i=1}^N p_i \log (p_i)}^{\text{Shannon entropy}} -
\overbrace{\beta \left(\sum_{i=1}^N E_i p_i - \langle E \rangle \right)}^
{\text{average energy constraint}} -
\overbrace{\mu \left(\sum_{i=1}^N  p_i - 1 \right)}^
{\text{normalization constraint}},
$$
where $N$ is the total number of microstates available to the system, and
$\beta$, and $\mu$ are the Lagrange multipliers associated with each of the
constraints. The next step consists on computing the gradient of this Lagrangian
which returns a vector of size $N$ where the $k^{\text{th}}$ entry is the
derivative of the Lagrangian with respect to $p_k$. But notice that all of these
derivatives will look the same. So taking one of these derivatives is enough.
We then take the derivative with respect to a particular $p_k$ and equate it
to zero, obtaining
$$
\frac{d\mathcal{L}}{d p_k} = -\log(p_k) - 1 - \lambda - \beta E_k = 0.
$$
Notice that all of the terms with $i\neq k$ disappear, leaving a simple
expression. Solving for $p_k$ gives
$$
p_k = \exp\left[1 - \lambda - E_k \right] = e^{1 - \lambda} e^{-\beta E_k}.
$$
Every single probability $p_k$ takes the same form. We substitute this
probability $p_k$ on our normalization constraint, obtaining
$$
\sum_{i=1}^N p_i = e^{1 - \lambda} \sum_{i=1}^N e^{-\beta E_i}=1.
$$
This tells us that the term $e^{1 - \lambda}$ is given by
$$
e^{1 - \lambda} = \frac{1}{\sum_{i=1}^Ne^{-\beta E_i}}.
$$
Therefore, the probability of microstate $i$ is given by
$$
P(E_i) = p_i = \frac{e^{-\beta E_i}}{\sum_{i=1}^N e^{-\beta E_i}},
$$
exactly the Boltzmann distribution. One can show why it is the case that our
Lagrange multiplier $\beta$ is exactly $1/k_BT$ as demanded by the thermodynamic
version of this distribution, but that is out of the scope for our purposes.
This section aims only to show the subtle and deep connection between
statistical mechanics and information theory. This connection suggests that part
of the unreasonable effectiveness of statistical mechanics might not come from
the physical basis of its core theory; but instead from the statistical
inference problem on which, given the limited information we have of any
thermodynamic system's microstate, entropy maximization gives us a recipe on
what the best guess for the probability distribution over the microstates is.

### Joint Uncertainty in an Uncertain World

Part of the complexity in understanding biological systems is that their
components form a network of interactions. This connectivity means that one part
of the organism's state depends on many other parts' states. For example, the
wild-type *lac* operon's expression depends on the conformation state of two
transcription factors: CRP and LacI. The state of these transcription factors
depends on the concentration of cyclic-AMP and allolactose, respectively. These
concentrations rely on the state of the environment and transporters'
availability to bring them into the cell. This chain of connections continues
indefinitely.

The mathematical language to express the dependence between two variables is
that of joint and conditional probability. Shannon's entropy (Eq.
$\ref{eq:shannon_result}$) can also be extended to account for dependence
between variables. To make the notation for this extension easier to follow, let
us use a different notation from now on. Let us express Shannon's entropy as
$$
H(m) = -\sum_m P(m) \log P(m),
\label{eq:shannon_x}
$$
where instead of giving a vector of probabilities $\mathbf{p}$ to the function
$H$, we now give it a random variable $m$. This notation is understood as: the
entropy is calculated over the distribution of possible values that $m$ can
take. If $m$ can take values $\{m_1, m_2, \ldots, m_n\}$, the probability of
obtaining $m = m_k$ is given by the function $P(m=m_k)$, which for brevity we
can write simply as $P(m_k)$. What Eq. $\ref{eq:shannon_x}$ is saying is: take
the random variable $m$ and all the possible values it can have; compute the
Shannon entropy by summing over the probability of all those values. In this
way, $H(m)$ is a shorthand for writing $H[P(m)]$.

With this notation in hand, let us think about two correlated random variables
$m$ and $p$. These could be the number of mRNAs and proteins in the cells, as
depicted in [@Fig:ch1_fig11](A). The *joint entropy* $H(m, p)$ measures the
uncertainty we have about the outcome of a pair of variables rather than a
single. All it takes is to sum over both variables on Eq. $\ref{eq:shannon_x}$
as
$$
H(m, p) = -\sum_m \sum_p P(m, p) \log P(m, p).
\label{eq:joint_entropy}
$$
Eq. $\ref{eq:joint_entropy}$ then does the same computation as Eq.
$\ref{eq:shannon_x}$, except that the sum is taken over all possible pairs of
random variables $m$ and $p$. But what if we get to observe the outcome of one
of the two variables (observing mRNA via RNA-seq, for example), can that tell us
something about the outcome of the other one? For this, we need to understand
the concept of conditional entropy.

![**Shannon's entropy for more than one random variable.** (A) Toy model of a
random process where mRNA (random variable $m$) is stochastically produced as a
Poisson process with a fixed mean. Proteins (random variable $p$) are also
stochastically produced as a Poisson process, but the mean depends on the number
of mRNAs. (B) Samples from the model presented in (A). The center plot shows the
joint distribution $P(m, p)$, while the edge histograms show the marginal
distributions $P(m)$ and $P(p)$. (C) Venn diagram of the relationship of
different information metrics. The [Python code
(`ch1_fig11.py`)](https://github.com/mrazomej/phd/blob/master/src/chapter_01/code/ch1_fig11.py)
used to generate this figure can be found on the thesis [GitHub
repository](https://github.com/mrazomej/phd).](ch1_fig11){#fig:ch1_fig11
short-caption="Shannon's entropy for more than one random variable"}
 
### Thinking Conditionally, a Condition for Thinking

In Joe Blitztein's excellent *Introduction to probability* [@Blitzstein2019], he
clarifies how conditional probability is one of the most powerful concepts in
probability theory. Through the concept of conditional probability, we can learn
whether or not two things are somehow correlated, allowing us from there to
dissect the nature of such correlation. Given the probabilistic nature of
Shannon's entropy, the power of conditional entropy is extended to the so-called
conditional entropy $H(p \mid m)$. Let us think of our two random variables $m$
and $p$ with a joint probability distribution $P(m, p)$. We can assume that the
outcome of both random variables is correlated for our mRNA-protein pair,
meaning that specific pairs of values are more likely to appear. If we observed
the outcome of one of the random variables and knew the correlation function
between random variables, our guess for the variable's value that we did not
observe would improve over a completely random choice. In our example, if we get
to observe that $m$ is a small (or large) number, we would suspect that $p$ is
also a small (or large) number, as shown in [@Fig:ch1_fig11](B). This means that
our uncertainty on the value of $p$ changed---it was reduced---upon observing
the value of $m$. The new uncertainty, i.e., the entropy of $p$ having learned
the value of $m$, averaged over all possible values of $m$, is computed as
$$
H(p \mid m) = - \sum_m \sum_p P(m) P(p \mid m) \log P(p \mid m),
\label{eq:conditional_entropy}
$$
where $P(p \mid m)$ is read as "probability of $p$ given that we observe $m$."
Finally, with all these concepts in hand, we can discuss the idea of information
in the Shannon sense.

### One Person's Entropy is Another Person's Information

So far, our discussion has focused on the concept of entropy. We first derived
the Shannon entropy from three basic principles that a metric of uncertainty
should satisfy. Then, we showed that one of the main statistical mechanics
results, i.e., the Boltzmann distribution, could be derived from maximizing this
entropy subject to certain constraints, suggesting that statistical mechanics
could be nothing more than an optimal statistical inference protocol, given
limited information. But no mention of information up to now. This intentional
omission is because we first needed to master the idea of entropy to understand
the mathematical definition of information. 

Recall that $H(p)$ quantifies the uncertainty about the outcome of the random
process that generates the value of the variable $p$. Furthermore, $H(p \mid m)$
quantifies the uncertainty about the outcome of the same variable, but this time
observing the outcome of the random variable $m$. In the worst-case scenario,
$m$ and $p$ are uncorrelated, and learning the value of $m$ does not tell us
anything about $p$. In that case, we then have that
$$
 H(p \mid m) = H(p)\;\; \text{ for $m$ and $p$ uncorrelated}.
$$
If $m$ and $p$ are correlated, as depicted in [@Fig:ch1_fig11](B), then the
uncertainty about $p$ is reduced upon learning the value of $m$, giving us a
general relationship between marginal and conditional entropy of the form
$$
H(p) \geq H(p \mid m).
$$
In this latter scenario, learning the value of $m$ reduced our uncertainty in
the possible value of $p$. This reduction in uncertainty agrees with an informal
definition of what "obtaining information" means. We can then define the mutual
information $I(m;p)$ between random variable $m$ and $p$ as the reduction in
uncertainty about the value of one of the random variables when we learn the
value of the other random variable. For our example in which we get to observe
the mRNA copy number, this would mean that the mutual information is computed as
$$
I(m; p) \equiv H(p) - H(p \mid m).
\label{eq:mutual_info_entropy}
$$
But the mutual information is symmetric, meaning that the information about the
outcome of one of the variables by observing the other variables is the same
when the roles of what we get to observe are inverted. This argument means that
we can mathematically show that
$$
I(m; p) = H(m) - H(m \mid p).
$$
This symmetry is why traditionally, the mutual information is written with a
semi-colon rather than a regular comma, indicating that the order of the
variables does not matter. To show the above symmetry, let us substitute the
definitions of the marginal conditional entropy. This substitution for Eq.
$\ref{eq:mutual_info_entropy}$ results in
$$
I(m; p) = 
- \sum_p P(p)\log P(p) - 
\left[ - \sum_m \sum_p P(m) P(p \mid m) \log P(p \mid m)\right].
\label{eq:mutual_info_probs}
$$
The trick is now to use the definition of conditional probability in the right
way. We know that the conditional probability is defined as
$$
P(p \mid m) \equiv \frac{P(m, p)}{P(m)}.
\label{eq:cond_prob}
$$
Furthermore, we know that we can obtain the probability $P(p)$ by marginalizing
the joint distribution $P(m, p)$ over all values of $m$. Mathematically this is
written as
$$
P(p) = \sum_m P(m, p).
\label{eq:marginalization}
$$
What Eq. $\ref{eq:marginalization}$ is stating is that to compute the
probability of observing value $p$ of our random variable, we can add the
probability of all pairs $m, p$ with the desired that have the desired value of
$p$. For Eq. $\ref{eq:mutual_info_probs}$, we substitute Eq.
$\ref{eq:marginalization}$ on the first term (outside of the $\log$) of the
right-hand side and Eq.~$\ref{eq:cond_prob}$ on the second term (in and outside
of the $\log$), obtaining
$$
I(m; p) = 
- \sum_p \left[ \sum_m P(m, p) \right]\log P(p)
+ \sum_m \sum_p P(m, p) \log \frac{P(m, p)}{P(m)}.
$$
Since the order of the sums do not matter, we can factorize the common terms
on the left-hand side and use the properties of logarithms to write
$$
I(m; p) = \sum_m \sum_p P(m, p) \log \frac{P(m, p)}{P(m)P(p)}.
$$
It is now easier to see that we would arrive at the same result if we started
with the opposite conditional entropy $P(m \mid p)$. These series of
manipulations where we write either joint or conditional entropies will become
handy in this thesis as we explore biophysical models of how to compute gene
expression input-output functions (more on that in Chapter 3).
[@Fig:ch1_fig11](C) shows a schematic representation of the relationship of all
the entropy-based quantities that we explored in this chapter. Although it is
impossible to cover an entire field in a short introduction, I hope this
intuitive explanation will suffice to understand the rest of the thesis.
