# pyBET_evo
PyBET_evo is a software implementation of the Bayesian Event Tree for Eruption Forecasting (BET_EF) (Marzocchi et al., 2008).  
Through this software the user can compute the long-term probability of unrest (Node 1), magmatic unrest (Node 2) and eruption (Node 3) when the volcano is in a quiet state, or the short-term probability of magmatic unrest and eruption when the volcano is undergoing unrest, by using the entropy-based method developed by Marzocchi et al. (2024).  
The peculiarity of this code is that the user can evaluate the evolution of the short-term probabilities on multiple time windows at the same time, which was not possible with the previous implementations (Tonini et al., 2015).  
An application of this code is shown in in Ferrara et al. (2025) to forecast the evolution of the recent phase of unrest at Campi Flegrei by defining anomalies thorugh experts' elicitation (Selva et al., 2012).  

# References
- Ferrara, S., Selva, J., Sandri, L., Marzocchi, & Elicitation VI Working Group (2025). Forecasting the evolution of the current unrest of Campi Flegrei by defining anomalies through experts’ elicitation. Annals of Geophysics, 68(1), V108-V108.  
- Marzocchi, W., Sandri, L., & Selva, J. (2008). BET_EF: a probabilistic tool for long-and short-term eruption forecasting. Bulletin of Volcanology, 70, 623-632.  
- Marzocchi, W., Sandri, L., Ferrara, S., & Selva, J. (2024). From the detection of monitoring anomalies to the probabilistic forecast of the evolution of volcanic unrest: an entropy-based approach. Bulletin of Volcanology, 86(1), 5.  
- Selva J., Marzocchi W., Papale P., Sandri L. (2012). Operational eruption forecasting at high-risk volcanoes: the case of Campi Flegrei, Naples. Journal of Applied Volcanology, 1, 1-14.
- Tonini, R., Sandri, L., & Thompson, M. A. (2015). PyBetVH: A Python tool for probabilistic volcanic hazard assessment and for generation of Bayesian hazard curves and maps. Computers & Geosciences, 79, 38-46.
