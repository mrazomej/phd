## Flow Cytometry {#sec:ch4_sec05}

In this section, we provide information regarding the equipment used to make
experimental measurements of the fold-change in gene expression in the interests
of transparency and reproducibility. We also provide a summary of our
unsupervised method of gating the flow cytometry measurements for consistency
between experimental runs.

### Equipment 

Due to past experience using the Miltenyi Biotec MACSQuant flow cytometer during
the Physiology summer course at the Marine Biological Laboratory, we used the
same flow cytometer for the formal measurements in this work graciously provided
by the Pamela Björkman lab at Caltech. All measurements were made using an
excitation wavelength of $488\,\text{nm}$ with an emission filter set of
525/$50\,\text{nm}$. This excitation wavelength provides approximately 40% of
the maximum YFP absorbance [@SpectraViewer], and this was found to be sufficient
for the purposes of these experiments. A useful feature of modern flow cytometry
is the high-sensitivity signal detection through the use of photomultiplier
tubes (PMT) whose response can be tuned by adjusting the voltage. Thus, the
voltage for the forward-scatter (FSC), side-scatter (SSC), and gene expression
measurements were tuned manually to maximize the dynamic range between
autofluorescence signal and maximal expression without losing the details of the
population distribution. Once these voltages were determined, they were used for
all subsequent measurements. Extremely low signal producing particles were
discarded before data storage by setting a basal voltage threshold, thus
removing the majority of spurious events. The various instrument settings for
data collection are given in [@tbl:ch4_tab01].

| **Laser**        | **Channel**                     | **Sensor Voltage** |
| ---------------- | ------------------------------- | ------------------ |
| $488\,\text{nm}$ | Forward-Scatter (FSC)           | $423\,\text{V}$    |
| $488\,\text{nm}$ | Side-Scatter (SSC)              | $537\,\text{V}$    |
| $488\,\text{nm}$ | Intensity (B1 Filter, 525/50nm) | $790\,\text{V}$    |
| $488\,\text{nm}$ | Trigger (debris threshold)      | $24.5\,\text{V}$   |
Table: **Instrument settings for data collection using the Miltenyi Biotec
MACSQuant flow cytometer.** All experimental measurements were collected using
these values. {#tbl:ch4_tab01}

### Experimental Measurement 

Prior to each day's experiments, the analyzer was calibrated using MACSQuant
Calibration Beads (Cat. No. 130-093-607) such that day-to-day experiments would
be comparable. A single data set consisted of seven bacterial strains, all
sharing the same operator, with varying repressor copy numbers ($R = 0$, 22, 60,
124, 260, 1220, and 1740), in addition to an autofluorescent strain, under
twelve IPTG concentrations. Data collection took place over two to three hours.
During this time, the cultures were held at approximately 4$^\circ$C by placing
the 96-well plate on a MACSQuant ice block. Because the ice block thawed over
the course of the experiment, the samples measured last were approximately at
room temperature. This means that samples may have grown slightly by the end of
the experiment. To confirm that this continued growth did not alter the measured
results, a subset of experiments were run in reverse meaning that the fully
induced cultures were measured first and the uninduced samples last. The plate
arrangements and corresponding fold-change measurements are shown in and ,
respectively. The measured fold-change values in the reverse ordered plate
appear to be drawn from the same distribution as those measured in the forward
order, meaning that any growth that might have taken place during the experiment
did not significantly affect the results. Both the forward and reverse data sets
were used in our analysis.

![**Plate arrangements for flow cytometry.** (A) Samples were measured primarily
in the forward arrangement with a subset of samples measured in reverse. The
black arrow indicates the order in which samples were processed by the flow
cytometer. (B) The experimentally measured fold-change values for the two sets
of plate arrangements show that samples measured in the forward arrangement
appear to be indistinguishable from those measured in reverse
order.](ch4_fig08){#fig:ch4_fig08 short-caption="Plate arrangements for flow
cytometry."}

### Unsupervised Gating 

Flow cytometry data will frequently include a number of spurious events or other
undesirable data points such as cell doublets and debris. The process of
restricting the collected data set to those data determined to be "real" is
commonly referred to as gating. These gates are typically drawn manually
[@Maecker2005] and restrict the data set to those points which display a high
degree of linear correlation between their forward-scatter (FSC) and
side-scatter (SSC). The development of unbiased and unsupervised methods of
drawing these gates is an active area of research [@Lo2008; @Aghaeepour2013].

For this study, we used an automatic unsupervised gating procedure to filter the
flow cytometry data based on the front and side-scattering values returned by
the MACSQuant flow cytometer. We assume that the region with highest density of
points in these two channels corresponds to single-cell measurements. Everything
extending outside of this region was discarded in order to exclude sources of
error such as cell clustering, particulates, or other spurious events.

In order to define the gated region we fit a two-dimensional Gaussian
function to the $\log_{10}$ forward-scattering (FSC) and the $\log_{10}$
side-scattering (SSC) data. We then kept a fraction $\alpha \in [0, 1]$
of the data by defining an elliptical region given by 
$$
\left(\boldsymbol{x} - \boldsymbol{\mu} \right)^T \boldsymbol{\Sigma}^{-1}
\left(\boldsymbol{x} - \boldsymbol{\mu} \right) \leq \chi^2_\alpha(p),
\label{eq:ch4_eq19}
$$
where $\boldsymbol{x}$ is the $2 \times 1$ vector containing the
$\log(\text{FSC})$ and $\log(\text{SSC})$, $\boldsymbol{\mu}$ is the $2 \times
1$ vector representing the mean values of $\log(\text{FSC})$ and
$\log(\text{SSC})$ as obtained from fitting a two-dimensional Gaussian to the
data, and $\boldsymbol{\Sigma}$ is the $2\times 2$ covariance matrix also
obtained from the Gaussian fit. $\chi^2_\alpha(p)$ is the quantile function for
probability $p$ of the chi-squared distribution with two degrees of freedom.
shows an example of different gating contours that would arise from different
values of $\alpha$ in . In this work, we chose $\alpha = 0.4$ which we deemed
was a sufficient constraint to minimize the noise in the data. As explained in
Appendix XXX we compared our high throughput flow cytometry data with single
cell microscopy, confirming that the automatic gating did not introduce
systematic biases to the analysis pipeline. The specific code where this gating
is implemented can be found in [GitHub
repository](https://github.com/RPGroup-PBoC/mwc_induction/blob/master/code/analysis/unsupervised_gating.ipynb).

![**Representative unsupervised gating contours.** Points indicate individual
flow cytometry measurements of forward scatter and side scatter. Colored points
indicate arbitrary gating contours ranging from 100\% ($\alpha = 1.0$) to 5\%
($\alpha = 0.05$). All measurements for this work were made computing the mean
fluorescence from the 40$^\text{th}$ percentile ($\alpha = 0.4$), shown as
orange points.](ch4_fig09){#fig:ch4_fig09 short-caption="Representative
unsupervised gating contours."}

### Comparison of Flow Cytometry with Other Methods 

Previous work from our lab experimentally determined fold-change for similar
simple repression constructs using a variety of different measurement methods
[@Garcia2011b; @Brewster2014]. Garcia and Phillips used the same background
strains as the ones used in this work, but gene expression was measured with
Miller assays based on colorimetric enzymatic reactions with the LacZ protein
[@Garcia2011c]. @Brewster2014 used a LacI dimer with the tetramerization region
replaced with an mCherry tag, where the fold-change was measured as the ratio of
the gene expression rate rather than a single snapshot of the gene output.

shows the comparison of these methods along with the flow cytometry method used
in this work. The consistency of these three readouts validates the quantitative
use of flow cytometry and unsupervised gating to determine the fold-change in
gene expression. However, one important caveat revealed by this figure is that
the sensitivity of flow cytometer measurements is not sufficient to accurately
determine the fold-change for the high repressor copy number strains in O1
without induction. Instead, a method with a large dynamic range such as the
Miller assay is needed to accurately resolve the fold-change at such low
expression levels.

![**Comparison of experimental methods to determine the fold-change.** The
fold-change in gene expression for equivalent simple-repression constructs has
been determined using three independent methods: flow cytometry (this work),
colorimetric miller assays [@Garcia2011c], and video microscopy [@Brewster2014].
All three methods give consistent results, although flow cytometry measurements
lose accuracy for fold-change less than $10^{-2}$. note that the repressor-dna
binding energies $\delta\varepsilon_{ra}$ used for the theoretical predictions
were determined in [@Garcia2011c].](ch4_fig10){#fig:ch4_fig10
short-caption="Comparison of experimental methods to determine the
fold-change."}

