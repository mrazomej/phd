## Sources of variability

In this section we will explore possible sources of variability that could
explain some of the discrepancies between our theoretical predictions and the
actual data. This section appeared as an intellectual exercise after our first
reviewing process. We thank the anonymous reviewer for bringing these points to
our attention as we keep on exploring the extent to which our general
theoretical framework can describe the data in terms of a small set of
parameters.

### Multiple gene copies

The first source of variability that we will explore is the fact that during the
cell cycle cells at the growth rate that we perform our experiments will have up
to two copies of the reporter gene . For this scenario we will explore two
possibilities:
1.  The case in which both promoters are independent of each other.
2.  The case in which the connection between both promoters due to a
    shared transcription factor pool is explicitly accounted for.

#### Independent promoters

We begin by defining $1 - f$ as the fraction of the time during the cell cycle
in which cells have only one copy of the reporter gene. This implies that $f$ is
the fraction of the time that cells spend with two copies of the gene. Let
$P_i(m)$ be the probability function of having $m$ mRNA molecules given that the
cell has $i \in \{1, 2 \}$ copies of the reporter gene. This means that $P(m)$,
the probability of any cell having $m$ mRNA molecules is given by 
$$
P(m) = (1 - f) P_1(m) + f P_2(m).
\label{eq:ch4_eq52}
$$

The average number of mRNA $\left\langle m \right\rangle$ is therefore given by 
$$
\left\langle m \right\rangle= \sum_{m=0}^\infty m P(m) =
(1 - f) \sum_{m=0}^\infty m P_1(m) +
f \sum_{m=0}^\infty m P_2(m).
\label{eq:ch4_eq53}
$$
This can be simply written as
$$
\left\langle m \right\rangle= (1 - f)
\left\langle m \right\rangle_1 + f \left\langle m \right\rangle_2
\label{eq:ch4_eq54}
$$
where $\left\langle m \right\rangle_i$ is the average number of mRNA molecules
for a cell with $i \in \{ 1 , 2 \}$ reporter gene copies.

It is perfectly reasonable to assume that $\left\langle m \right\rangle_2 = 2
\left\langle m \right\rangle_1$ if the promoter does not subtract much of the
cellular resources. This would imply that 
$$
\begin{split}
\left\langle m \right\rangle &= (1 - f) 
\left\langle m \right\rangle_1 + 2 f \left\langle m \right\rangle_1,\\
&= (1 + f) \left\langle m \right\rangle_1.
\end{split}
\label{eq:ch4_eq55}
$$

From a kinetic model of this transcriptional regulation we know that the mean
expression level for a single unregulated promoter is given by [@Jones2014a]
$$
\left\langle m \right\rangle_1(R = 0) = \frac{r}{\gamma},
\label{eq:ch4_eq56}
$$
where $r$ is the mRNA production rate and $\gamma$ is the mRNA degradation rate.
For a single regulated promoter then the mean mRNA copy number is given by
$$
\begin{split}
\left\langle m \right\rangle_1 (R > 0) &= \frac{r}{\gamma} \times 
\text{fold-change},\\
& = \frac{r}{\gamma} \times \frac{1}{
1 + \frac{R p_{\text{act}}}{N_{NS}} e^{-\beta \Delta \varepsilon_{RA}}}.
\end{split}
\label{eq:ch4_eq57}
$$

This means that when we compute the fold-change for the case of
independent promoters we have
$$
\begin{split}
\text{fold-change} &= \frac{\left\langle m \right\rangle(R > 0)}{
\left\langle m \right\rangle(R = 0)}, \\
&= \frac{(1 + f) \left\langle m \right\rangle_1 (R > 0)}{
(1 + f) \left\langle m \right\rangle_1 (R = 0)}.
\end{split}
\label{eq:ch4_eq58}
$$
Since the $(1 + f)$ factor cancels, the expression for fold-change is not
altered obtaining
$$
\text{fold-change} = \frac{1}{1 + \frac{R p_{\text{act}}}{N_{NS}}
e^{-\beta \Delta\varepsilon_{RA}}}.
\label{eq:ch4_eq59}
$$

#### Dependent promoters

When we relax the promoters independence assumption we can use the grand
canonical ensemble formulation as described in [@Weinert2014]. In this
description the fold-change equation is given by
$$
\text{fold-change} = \frac{1}{1 + \lambda_r e^{-\beta\Delta\varepsilon_{RA}}},
\label{eq:ch4_eq60}
$$
where $\lambda_r$ is the fugacity of the repressor. The value of this fugacity
is obtained by taking into account all of the repressor reservoirs considered in
the system. In our case there are two repressor reservoirs: repressors bound to
specific binding sites $\left\langle R_S \right\rangle$ and repressors bound to
non-specific binding sites $\left\langle R_{NS} \right\rangle$. These two
reservoirs are connected through the constraint
$$
\left\langle R_{\text{tot}} \right\rangle = 
\left\langle R_S \right\rangle + \left\langle R_{NS} \right\rangle,
\label{eq:ch4_eq61}
$$
where $\left\langle R_{\text{tot}} \right\rangle$ is the repressor copy number
as determined by an independent method -- quantitative westerns as in
[@Garcia2011c] or binomial partitioning as in [@Brewster2014] --. As shown in
[@Weinert2014] the repressor reservoir expressions are given by
$$
\left\langle R_S \right\rangle= N_{S}
\frac{\lambda_r e{-\beta\Delta\varepsilon_{RA}}}{
1 + \lambda_r e^{-\beta\Delta\varepsilon_{RA}}},
\label{eq:ch4_eq62}
$$
and
$$
\left\langle R_{NS} \right\rangle = N_{NS}\frac{\lambda_r}{1 + \lambda_r},
\label{eq:ch4_eq63}
$$
where $N_{S}$ is the number of specific binding sites for the repressor, i.e.
the number of promoters, and $N_{NS}$ is the number of non-specific binding
sites.

Since $\left\langle R_{\text{tot}} \right\rangle$ is measured using independent
methods we can constrain the value of $\lambda_r$ using in combination with and
XX. This is
$$
\left\langle R_{tot} \right\rangle = 
N_{S} \frac{ \lambda_r \exp( -\beta \Delta \varepsilon_{RA} ) }
{1 + \lambda_r \exp( -\beta \Delta \varepsilon_{RA} ) } 
+
N_{NS}\frac{\lambda_r}{1 + \lambda_r},
\label{eq:ch4_eq64}
$$
<!-- Distributing terms gives a second order polynomial on $\lambda_r$ of the form
$$
\lambda_r^2 e^{-\beta\Delta\varepsilon_{RA}} (\left\langle R_{\text{tot}}
\right\rangle- N_{S}- N_{NS}) +
\lambda_r (\left\langle R_{\text{tot}} \right\rangle +
\left\langle R_{\text{tot}} \right\rangle
e^{-\beta\Delta\varepsilon_{RA}} - N_{S}
e^{-\beta\Delta\varepsilon_{RA}} - N_{NS}) +
\left\langle R_{\text{tot}} \right\rangle= 0.
\label{eq:ch4_eq65}
$$

For the quadratic term in $\lambda_r$ we note that $N_{NS}\gg \left\langle
R_{\text{tot}} \right\rangle, N_{S}$, therefore we can approximate it as
$$
\lambda_r^2 e^{-\beta\Delta\varepsilon_{RA}}
(\left\langle R_{\text{tot}} \right\rangle- N_{S}- N_{NS}) \approx
\lambda_r^2 e^{-\beta\Delta\varepsilon_{RA}} N_{NS}.
\label{eq:ch4_eq66}
$$

<!-- For the linear term on $\lambda_r$ only the terms $\left\langle R_{\text{tot}}
\right\rangle e^{-\beta\Delta\varepsilon_{RA}})$ and
$N_{S}e^{-\beta\Delta\varepsilon_{RA}}$ are of the same order of magnitude as
$N_{NS}$, so we have
$$
\lambda_r (\left\langle R_{\text{tot}} \right\rangle +
\left\langle R_{\text{tot}} \right\rangle 
e^{-\beta\Delta\varepsilon_{RA}} - N_{S}e^{-\beta\Delta\varepsilon_{RA}} - 
N_{NS}) 
\approx \lambda_r \left[ e^{-\beta\Delta\varepsilon_{RA}} 
\left(\left\langle R_{\text{tot}} \right\rangle- N_{S}\right)  - N_{NS}\right].
\label{eq:ch4_eq67}
$$
With these simplifications in hand is a matter of finding the positive root of
this second order polynomial to obtain the value for the fugacity. -->

<!-- Using this formulation the mean mRNA copy number is then given by
$$
\left\langle m \right\rangle = (1 - f)
\left\langle m \right\rangle_1 + f \left\langle m \right\rangle_2.
\label{eq:ch4_eq68}
$$

If we again assume that having two copies of the promoter produces twice as much
as one promoter this is given by
$$
\left\langle m \right\rangle= (1 - f) \frac{r}{\gamma} \times
\text{fold-change}(N_{S}= 1) +
f \frac{(2r)}{\gamma} \times \text{fold-change}(N_{S}= 2),
\label{eq:ch4_eq69}
$$
where we explicitly state that the fold-change function depends on the number of
specific binding sites. -->

<!-- For the fold-change equation of the total mean number of mRNA we then have 
$$
\text{fold-change} =
\frac{ \frac{r}{\gamma} \left[ (1 - f) \times \text{fold-change}(Ns = 1) +
2 f \times \text{fold-change}(N_{S}=2) \right]
}{
\frac{r}{\gamma} \left[ (1 - f) + 2 f \right]. }
\label{eq:ch4_eq70}
$$
When substituting the definition of fold-change we obtain our
final expression 
$$
\text{fold-change} = \frac{1}{(1 + f)} \left[
(1 - f) \frac{1}{1 + \frac{R}{N_{NS}} e^{-\beta\Delta\varepsilon_{RA}}} +
2 f \frac{1}{1 + \lambda_r e^{-\beta\Delta\varepsilon_{RA}}} \right].
\label{eq:ch4_eq71}
$$

Using we can evaluate the effect of having two promoters for a fraction of the
cell cycle. For this we consider that the non-specific number of binding sites
$N_{NS}$ doubles for the case in which there are two copies of the promoter
since the reporter construct is close to the end of the replication fork. shows
the comparison of evaluating the single promoter fold-change equation (solid
lines) and (dotted lines). We can see that for the strongest binding sites O1
(panel A) and Oid (panel D) the prediction for the high inducer concentration
are significantly worse for the two promoter model model than for the single
promoter case. For the intermediate binding site level (O2 panel B) the
difference are not that significant compared with the single promoter case.
Finally for the case of the weakest binding site we see that this assumption
alleviates some of the discrepancies between theory and prediction; nevertheless
one could argue that similar improvements can be accounted for when performing a
global fit of all parameters as described in appendix
[\[appendix_global_fit\]](#appendix_global_fit){reference-type="ref"
reference="appendix_global_fit"}.

![**Comparison of the effect of a single vs a double promoter model.** Four
different binding sites for the repressor. Solid lines shows the predictions
made with the single-promoter model, dotted lines show the predictions of the
model that includes a fraction $f = 1/3$ of the cell cycle in which the cell has
two copies of the promoter.](ch4_fig30){#fig:ch4_fig30 short-caption="Comparison
of the effect of a single vs a double promoter model"}
