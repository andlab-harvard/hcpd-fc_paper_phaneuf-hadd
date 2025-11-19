# hcpd-fc_paper_phaneuf-hadd

Data and Code

*to accompany*

Phaneuf-Hadd, C.V., Heffer, T., Mair, P., & Somerville, L.H. (under review). Charting Age-Related Change in the Architecture of Fluid Cognition. *Child Development*. doi: PENDING.

## Developer Contact Information

Github profile: https://github.com/cphaneuf

Email: (current) cphaneuf@g.harvard.edu, (permanent) cphaneuf@umich.edu

## Contents

### data/ directory

Contains demographic and NIH Toolbox Cognition Battery data from the cross-sectional arm of the Human Connectome Project in Development (HCPD).

### analyses/ directory

*utilities.R* defines variables and functions to be shared across scripts.

*demog.R* takes demographic data inputs from data/ and writes outputs to results/demog/.

*gam.R* takes NIH Toolbox Cognition Battery data inputs from data/ and writes outputs (from generalized additive models) to results/gam/.

*psychom.R* takes NIH Toolbox Cognition Battery data inputs from data/ and writes outputs (from multidimensional scaling, factor analyses, correlation matrices, and correlation networks) to results/psychom/.

### results/ directory

Contains text and png file outputs from scripts in analyses/, sorted by analysis type.

### annotated_figs/ directory

Contains annotated_figs.pptx, which annotates several figures beyond the limits of R.
