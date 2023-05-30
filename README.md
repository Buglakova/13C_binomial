# 13C_binomial
Fit binomial model to the isotopologue distribution for fatty acids as described in the manuscript  

This analysis doesn't require any unusual libraries, still `environment.yml` to keep track of Python and packages versions.  

Function `fit_binomial` in `fit_binomial.py` takes as input numpy array of normalized and corrected for natural isotope abundance intensities for each isotopologue peak of the fatty acid and the model type. Model can be either `C16` or `C18`. `C18` takes into account that M+2 peak can come from the unlabeled palmitate + $2^{13}C$.

`example_isocorrector_analysis.ipynb` shows an example of how to analyze data coming from Isocorrector software.

## Output explanation
Example output for data `data/IsoCorrectoR_ACLY.csv` is in `data/acly_results`

- Plots show original data as bars and prediction of the fitted model as a line

- `fit_binomial_prediction.csv` contains corresponding predictions for each sample (table of the same format as input but with different values)

- `fit_binomial_result.csv` contains result of fit for each sample  
  model: which model was used, C16 or C18  
  fit_success: 1 if optimization procedure reached convergence condition  
  uptake: model parameter uptake  
  uptake_palmitate: model parameter palmitate uptake, only for C18  
  p: parameter of binomial distribution  
  mean: mean of isotopologue distribution ($sum(i * I(M + i))/sum(I(M+i))$), with M+0 set to 0 and additionally M+2 set to 0 for C18  
  p_mean: parameter of binomial distribution estimated from the distribution mean (for binomial $mean = pn$)
