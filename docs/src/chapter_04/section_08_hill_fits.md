## Alternate Characterizations of Induction 

In this section we discuss a different way to describe the induction data,
namely, through using the conventional Hill approach. We first demonstrate how
using a Hill function to characterize a single induction curve enables us to
extract features (such as the midpoint and sharpness) of that single response,
but precludes any predictions of the other seventeen strains. We then discuss
how a thermodynamic model of simple repression coupled with a Hill approach to
the induction response can both characterize an induction profile and predict
the response of all eighteen strains, although we argue that such a description
provides no insight into the allosteric nature of the protein and how mutations
to the repressor would affect induction. We conclude the section by discussing
the differences between such a model and the statistical mechanical model used
in the main text.

### Fitting Induction Curves using a Hill Function Approach

The Hill equation is a phenomenological function commonly used to describe data
with a sigmoidal profile [@Murphy2007; @Murphy2010; @Rogers2015]. Its simplicity
and ability to estimate the cooperativity of a system (through the Hill
coefficient) has led to its widespread use in many domains of biology
[@Frank2013]. Nevertheless, the Hill function is often criticized as a
physically unrealistic model and the extracted Hill coefficient is often
difficult to contextualize in the physics of a system [@Weiss1997]. In the
present work, we note that a Hill function, even if it is only used because of
its simplicity, presents no mechanism to understand how a regulatory system's
behavior will change if physical parameters such as repressor copy number or
operator binding energy are varied. In addition, the Hill equation provides no
foundation to explore how mutating the repressor (e.g., at its inducer-binding
interface) would modify its induction profile, although statistical mechanical
models have proved capable of characterizing such scenarios [@Keymer2006;
@Swem2008; @Einav2016].

Consider the general Hill equation for a single induction profile given
by 
$$
\text{fold-change} = (\text{leakiness}) + (\text{dynamic range}) 
\frac{\left( \frac{c}{K} \right)^n}{1 + \left( \frac{c}{K} \right)^n},
$${#eq:ch4_eq27}
where, as in the main text, the leakiness represents the minimum fold-change,
the dynamic range represents the difference between the maximum and minimum
fold-change, $K$ is the repressor-inducer dissociation constant, and $n$ denotes
the Hill coefficient that characterizes the sharpness of the curve ($n > 1$
signifies positive cooperativity, $n = 1$ denotes no cooperativity, and $n < 1$
represents negative cooperativity). shows how the individual induction profiles
can be fit (using the same Bayesian methods as described in Section "") to this
Hill response, yielding a similar response to that shown in . However,
characterizing the induction response in this manner is unsatisfactory because
each curve must be fit independently thus removing our predictive power for
other repressor copy numbers and binding sites.

The fitted parameters obtained from this approach are shown in . These are
rather unsatisfactory because they do not clearly reflect the properties of the
physical system under consideration. For example, the dissociation constant $K$
between LacI and inducer should not be affected by either the copy number of the
repressor or the DNA binding energy, and yet we see upward trends as $R$ is
increased or the binding energy is decreased. Here, the $K$ parameter ultimately
describes the midpoint of the induction curve and therefore cannot strictly be
considered a dissociation constant. Similarly, the Hill coefficient $n$ does not
directly represent the cooperativity between the repressor and the inducer as
the molecular details of the copy number and DNA binding strength are subsumed
in this parameter as well. While the leakiness and dynamic range describe
important phenotypic properties of the induction response, this Hill approach
leaves us with no means to predict them for other strains. In summary, the Hill
equation cannot predict how an induction profile varies with repressor copy
number, operator binding energy, or how mutations will alter the induction
profile. To that end, we turn to a more sophisticated approach where we use the
Hill function to describe the available fraction of repressor as a function of
inducer concentration.

![**Hill function and MWC analysis of each induction profile.** Data for each
individual strain was fit to the general Hill function in XXX. (A) strains with
O1 binding site, (B) strains with O2 binding site, and (C) strains with O3
binding site. Shaded regions indicate the bounds of the 95\% credible
region.](ch4_fig16){#fig:ch4_fig16 short-caption="Hill function and MWC analysis
of each induction profile"}

![**Parameter values for the Hill equation fit to each individual titration.**
The resulting fit parameters from the Hill function fits of [@Fig:ch4_fig16] are
summarized. The large parameter intervals for many of the O3 strains are due to
the flatter induction profile (as seen by its smaller dynamic range), and the
ability for a large range of $K$ and $n$ values to describe the
data.](ch4_fig17){#fig:ch4_fig17 short-caption="Parameter values for the Hill
equation fit to each individual titration"}

### Fitting Induction Curves using a Combination Thermodynamic Model and Hill Function Approach

Motivated by the inability in the previous section to characterize all eighteen
strains using the Hill function with a single set of parameters, here we combine
the Hill approach with a thermodynamic model of simple repression to garner
predictive power. More specifically, we will use the thermodynamic model in but
substitute the statistical model in with the phenomenological Hill function XXX.

Following , fold-change is given by
$$
\text{fold-change} = \left( 1 + p_A(c) \frac{R}{N_{NS}}e^{-\beta
\Delta\varepsilon_{RA}} \right)^{-1}
$${#eq:ch4_eq28}
where the Hill function
$$
p_A(c) = p_A^{\text{max}} - p_A^{\text{range}}
\frac{\left( \frac{c}{K_D} \right)^n}{1 + \left( \frac{c}{K_D} \right)^n}
$${#eq:ch4_eq29}
represents the fraction of repressors in the allosterically active state, with
$p_A^{\text{max}}$ denoting the fraction of active repressors in the absence of
inducer and $p_A^{\text{max}} - p_A^{\text{range}}$ the minimum fraction of
active repressors in the presence of saturating inducer. The Hill function
characterizes the inducer-repressor binding while the thermodynamic model with
the known constants $R$, $N_{NS}$, and $\Delta\varepsilon_{RA}$ describes how
the induction profile changes with repressor copy number and repressor-operator
binding energy.

As in the main text, we can fit the four Hill parameters -- the vertical shift
and stretch parameters $p_A^{\text{max}}$ and $p_A^{\text{range}}$, the Hill
coefficient $n$, and the inducer-repressor dissociation constant $K_D$--for a
single induction curve and then use the fully characterized to describe the
response of each of the eighteen strains. shows this process carried out by
fitting the O2 $R=260$ strain (white circles in Panel XXX) and predicting the
behavior of the remaining seventeen strains.

![**A thermodynamic model coupled with a Hill analysis can characterize
induction.** Combining a thermodynamic model of simple repression with the Hill
function to characterize the repressor-inducer binding successfully
characterizes the induction profiles of all eighteen strains. As in the main
text, data was only fit for the O2 $R=260$ strain using and the parameters
$p_A^{\text{max}} =0.90^{+0.03}_{-0.01}$, $p_A^{\text{range}} =
-0.90^{+0.02}_{-0.03}$, $n = 1.6_{-0.1}^{+0.2}$, and $K_D = 4^{+2}_{-1} \times
10^{-6}\,\text{M}$. Shaded regions indicate bounds of the 95% credible
region.](ch4_fig18){#fig:ch4_fig18 short-caption="A thermodynamic model coupled
with a Hill analysis can characterize induction"}

Although the curves in are nearly identical to those in (which were made using
the MWC model ), we stress that the Hill function approach is more complex than
the MWC model (containing four parameters instead of three) and it obscures the
relationships to the physical parameters of the system. For example, it is not
clear whether the fit parameter $K_D = 4^{+2}_{-1} \times 10^{-6}\,\text{M}$
relays the dissociation constant between the inducer and active-state repressor,
between the inducer and the inactive-state repressor, or some mix of the two
quantities.

In addition, the MWC model naturally suggests further quantitative tests for the
fold-change relationship. For example, mutating the repressor's inducer binding
site would likely alter the repressor-inducer dissociation constants $K_A$ and
$K_I$, and it would be interesting to find out if such mutations also modify the
allosteric energy difference $\Delta\varepsilon_{AI}$ between the repressor's
active and inactive conformations. For our purposes, the Hill function falls
short of the connection to the physics of the system and provides no intuition
about how transcription depends upon such mutations. For these reasons, we
present the thermodynamic model coupled with the statistical mechanical MWC
model approach in the paper.
