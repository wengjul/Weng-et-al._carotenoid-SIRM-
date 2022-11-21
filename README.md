# Weng-et-al._carotenoid-SIRM-
Single Cell Raman Microspectroscopy Evaluation R scripts

These codes are designed to evaluate single cell Raman microspectroscopy data.

All files contain visualization sections.

Weng-et-al._SCRM_general is a general evaluation script with steps like spectral cropping, baseline correction, min-max scaling, averaging etc.
Input are individual single cell Raman spectra. A Si waver spectrum is used for wavenumber correction.
Output are average spectra for each data folder.


Weng-et-al._mean spectra overview_long and Weng-et-al._mean spectra overview_short are examples for scripts for visualization of stacked average spectra for biomass SCRM and carotenoid SCRR spectra, respectively.

Weng-et-al._PCA_short is an example for a principal component analysis script.

Weng-et-al._signal position evaluation_gaussian is an example for a Gaussian fit based algorithm for determination of signal positions

Weng-et-al._signal integral evaluation_gaussian is an example for a Gaussian fit based algorithm for determination of signal integrals and integral ratios

Further information and citation of the used R packages can be found in the supporting information of Weng et al.


R version 4.2.1 (2022-06-23 ucrt)
