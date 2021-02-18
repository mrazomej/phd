## Definition of the non-specific background $N_{NS}$

In this section we will explore the definition of the non-specific background
$N_{NS}$. As raised by an anonymous reviewer the nature of this parameter seems
to raise some controversy on what the right value should be, or whether or not
the arbitrary definition of its value should also be applied to the
$\Delta\varepsilon_{AI}$ parameter.

Specifically during the first round a reviewer did not like the idea that the
value of $N_{NS} = 4.6 \times 10^6$ assumed that the entirety of the genome was
available for non-specific binding of the repressor. We will consider how
reasonable this is at the end of the section. However, As we will show first,
the specific value of $N_{NS}$ is analogous to the zero potential energy or the
reference concentration state. It is only the free energy differences that
matter at the end of the day. For the second round of reviews the same reviewer
was willing to agree on our point if and only if we were to acknowledge that
other parameters such as $\Delta\varepsilon_{AI}$, the free energy difference
between the active and inactive state of the repressor, also had an arbitrary
definition that could be set to any value. In this section we will show as well
that such a statement is an erroneous interpretation of the parameters. This
free energy difference value cannot be re-defined to take any value if one is to
be consistent with the experimental data.

Let us start by showing why the specific value of $N_{NS}$ is not the key
variable. Under the weak promoter approximation the fold-change equation is
equivalent to a two-state Fermi function of having the promoter being occupied
by a repressor or having an empty promoter. This is
$$
\text{fold-change} \rightarrow p^r_{bound} = \frac{1}{ 1 +
\frac{R}{N_{NS}} e^{-\beta\Delta\varepsilon_{RA}}}.
\label{eq:ch4_eq41}
$$
This expression can be rewritten as
$$
p^r_{bound} = \frac{1}{ 1 + e^{-\beta\Delta E}},
\label{eq:ch4_eq42}
$$
where $\Delta E$ is the free energy difference between the empty and occupied
promoter. This definition implies that
$$
\Delta E\equiv \overbrace{\Delta\varepsilon_{RA}}^{\text{enthalpic term}}
  - \overbrace{k_BT \ln \left( \frac{R}{N_{NS}} \right)}^{\text{entropic term}}.
\label{eq:ch4_eq43}
$$

Given that the parameter $\Delta\varepsilon_{RA}$ is inferred rather than
directly measured, this puts us in the position of being able to redefine
$N_{NS}$ at will as long as $\Delta E$ is in accordance to the experimental
data. In other words the parameter that matters is the free energy difference,
rather than its individual components. For example, if for a given operator and
a given repressor copy number we choose a different value of $N_{NS}$, it still
should hold true that
$$
\Delta E= \Delta\varepsilon_{RA}' - k_BT \ln \left( \frac{R}{N_{NS}'} \right),
\label{eq:ch4_eq44}
$$
where $N_{NS}'$ is the changed value of the non-specific background and
$\Delta\varepsilon_{RA}'$ is a different value for the repressor binding energy
that compensates for the difference in the non-specific background.

Let $N_{NS}' \equiv \alpha N_{NS}$, since the value of $\Delta E$ has to be
preserved it should be true that
$$
\Delta E= \Delta\varepsilon_{RA}' - 
k_BT \ln \left( \frac{R}{\alpha N_{NS}} \right)
= \Delta\varepsilon_{RA} - k_BT \ln \left( \frac{R}{N_{NS}} \right).
\label{eq:ch4_eq45}
$$
Solving for $\Delta\varepsilon_{RA}'$ gives 
$$
\begin{split}
\Delta\varepsilon_{RA}' &= \Delta\varepsilon_{RA} + 
k_BT \ln \left( \frac{N_{NS}}{\alpha N_{NS}} \right)\\
&= \Delta\varepsilon_{RA} - k_BT \ln \alpha.
\end{split}
\label{eq:ch4_eq46}
$$

Eq. $\ref{eq:ch4_eRA_redef}$ implies that we can redefine $N_{NS}$ to be any value as long as
$\Delta\varepsilon_{RA}$ compensates to maintain the value of $\Delta E$. This
statement holds true whether we are considering a single promoter or multiple
promoters. The same cannot be said about the $\Delta\varepsilon_{AI}$ parameter.
The parameter $\Delta\varepsilon_{AI}$ by itself sets the fraction of inactive
repressors in the absence of inducer via
$$
p_{act} = \frac{1}{1 + e^{-\beta\Delta\varepsilon_{AI}}},
\label{eq:ch4_eq47}
$$
where we have again a Fermi function for a two-state system in which the
repressor can be in an active or inactive state.

As shown before, the reason why we could define $N_{NS}$ to be any value is
because the parameter that matters is itself $\Delta E$ the free energy
difference. Therefore the repressor binding energy $\Delta\varepsilon_{RA}$
could compensate for changes in the value of $N_{NS}$. For the case of
$\Delta\varepsilon_{AI}$ tells us that $\Delta\varepsilon_{AI}$ has no entropic
term that can be compensated with an enthalpic term, or vice versa.

One could argue that for the case of a single promoter the fold-change equation
does allow this parameter to be redefined in an arbitrary way since the full
equation in the absence of inducer can be written as
$$
\text{fold-change} = \frac{1}{
1 + \left( \frac{1}{1 + e^{-\beta\Delta\varepsilon_{AI}}} \right)
\frac{R}{N_{NS}} e^{-\beta\Delta\varepsilon_{RA}}}.
\label{eq:ch4_eq48}
$$
So when we define the free energy $\Delta E$ we would include an extra term of
the form
$$
\Delta E= \Delta\varepsilon_{RA} - 
k_BT \left[ \ln \left( \frac{R}{N_{NS}} \right) +
\ln \left( \frac{1}{1 + e^{-\beta\Delta\varepsilon_{AI}}} \right) \right].
\label{eq:ch4_eq49}
$$
If we were only to use the statement brought up by the anonymous reviewer would
be true since changes in $\Delta\varepsilon_{AI}$ could be compensated by
changes in $\Delta\varepsilon_{RA}$ or $N_{NS}$. But as specified in appendix
[\[AppendixModel\]](#AppendixModel){reference-type="ref"
reference="AppendixModel"} this is not the case for cells with multiple
promoters.

The case of multiple promoters can be handled using the Canonical ensemble as in
or using the Grand Canonical ensemble as detailed in [@Weinert2014]. Our point
is more clearly seen in the case of the Canonical ensemble. Under this formalism
the fold-change equation is given by [@Brewster2014]
$$
\text{fold-change} = \frac{
\sum_{m=0}^{\min (N,R)} \frac{R!}{(N_{NS})^m (R - m)!}
{N \choose m} e^{-\beta m \Delta\varepsilon_{RA}}(N - m)
}{
N \sum_{m=0}^{\min (N,R)} \frac{R!}{(N_{NS})^m (R - m)!}
{N \choose m} e^{-\beta m \Delta\varepsilon_{RA}}
},
\label{eq:ch4_eq50}
$$
where $N$ is the number of promoters. Notice that we can group the terms
including $N_{NS}$ and $\Delta\varepsilon_{RA}$ as 
$$
\text{fold-change} = \frac{
\sum_{m=0}^{\min (N,R)} \frac{R!}{(R - m)!}
{N \choose m} \left(
\frac{e^{-\beta \Delta\varepsilon_{RA}}}{N_{NS}}\right)^m (N - m)
}{
N \sum_{m=0}^{\min (N,R)} \frac{R!}{(R - m)!}
{N \choose m} \left(\frac{e^{-\beta \Delta\varepsilon_{RA}}}{N_{NS}}\right)^m},
\label{eq:ch4_eq51}
$$
to highlight that it is a combination of these two parameters that matters,
rather than their individual values. For the case of the
$\Delta\varepsilon_{AI}$ parameter this is not the case. Every term containing
$R$ on is effectively multiplied by . Since these terms are included inside the
factorials it is not true that a simple compensation by the other parameters
allow us to define $\Delta\varepsilon_{AI}$ to be any value. Therefore as
defined in appendix [\[AppendixModel\]](#AppendixModel){reference-type="ref"
reference="AppendixModel"} the parameter $\Delta\varepsilon_{AI}$ can be
independently inferred using multiple promoter measurements of fold change.

As a final note, we can also check whether $N_{NS} = 4.6 \times 10^6$ is at all
a reasonable value to use. One potential point of concern is whether the
chromosomal DNA is occupied by other transcription factors that may reduce the
availability of the DNA for repressor or RNAP to bind. Here we consider data
from a recent census of protein abundance across the *E. coli* genome. In that
work, Schmidt *et al.* [@Schmidt2015] measured the protein copy number across
more than half the coding genes (greater than 95% by total protein mass). During
exponential growth in M9 minimal media with 0.5 % glucose, they find that about
6 % of the protein mass, or 311,000 monomer copies per cell, are proteins such
as transcription factors that will be bound to the DNA (about two-thirds of
these are nucleoid-associated proteins such as HNS and HU).

To make a simple estimate of DNA occupancy, let us assume that all transcription
factors bind DNA as dimers and occupy a DNA length of 15 bp (this appears to
vary from 7 bp to 38 bp in *E. coli* on RegulonDB [@Gama-Castro2016]), we find
that about 2.3 kbp or about half of the genome will be occupied. In the most
extreme case we could assume that this fraction is totally inaccessible, which
would reduce $N_NS$ by a factor of about $2$. Applying this to equation
[\[eq_eRA_redef\]](#eq_eRA_redef){reference-type="ref"
reference="eq_eRA_redef"}, we see that this has a negligible effect on the
actual binding energy that we would infer, and only corresponds to a change in
energy $\varepsilon_{RA}$ by about 0.7 $k_B T$.
