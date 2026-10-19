# spooky_correlation
This is a model of comprehensively tallying possible pairwise correlations of two quasars such that they could be exploited to mimic a violation of a Bell-inequality test. It reports the calculations and reproduces the outputs and plots of the closed-form versions in the journal paper. These are, effectively, a prediction of finding a pair of quasars that are correlated: one brightens, another dims - or vice-versa.  Various flavours of those calculations can be tried, to verify that they are mathematically and statistically correct as they appear in the paper.

Eric Steinbring, 2 January 2020; revised 21 September 2026

To run the code you need IDL and the astrolib libraries.  Download the PRO file and start IDL. The value of M (and R, if forcing a particular "residual" correlation) can be entered into the parameter-table, which is in the preamble of the code. The defaults are M=0.87 and R=0.67, for the holographic quantum-foam inspired result. Then, on the command line, type:

IDL>.compile spooky_correlation

IDL>spooky_correlation

And hit "enter."  It will generate the simulated distribution of flux-difference probabilities for given residual values M and R, and the plots of those as they appear in the paper.

You are free to use and modify this code as you wish.  If you do so, please cite:

Steinbring, E., "Detectability of Spooky Quasar-Flux Correlation on the Sky", Galaxies, 2026
