# Tutorial mmgsem
Repository containing tutorial scripts for the `mmgsem` R package. 

## Overview tutorial
The R script "overview_tutorial" is the main file accompanying the tutorial paper *"Using Mixture Multigroup Structural Equation Modeling to Compare Structural Relations Across Many Groups: An Introduction to the mmgsem R Package"*. It contains a step-by-step tutorial from measurement model estimation to hypothesis testing.

## Notes
### Appendix B Syntax Mistake
In *Appendix B: Estimating a Categorical Confirmatory Factor Analysis in `lavaan`* of the tutorial paper, there is a minor mistake when estimating the final model using `mmgsem()`. We skipped the `ordered` argument, which must be set to `TRUE`. The corrected snippet of that Appendix should look like this:

```r
mmgsem_fit <- mmgsem::mmgsem(
  dat     = data,
  S1      = S1_model,
  S2      = S2_model,
  group   = "group",
  nclus   = 4,
  seed    = 1,
  s1_type = "lavaan",
  s1_fit  = MM_fit_ordinal,
  ordered = TRUE
)
```

### Change Log
2026-10-28. The compute_se() function's default behavior changed. Specifically, the tolerance level for running ginv() on the Hessian matrix was set to 1e-05 for stability. This change resulted in slightly different pre-computed standard errors (i.e., `standard_errors` in the tutorial), which do not match the tutorial's results. For reproducibility, the deprecated object is now stored as `standard_errors_deprecated`, and can be accessed using `standard_errors <- mmgsem::standard_errors_deprecated`.

# References
Perez Alonso, A. F., Vermunt, J. K., Rosseel, Y., Zhao, H., & De Roover, K. (2026). Using Mixture Multigroup Structural Equation Modeling to Compare Structural Relations Across Many Groups: An Introduction to the mmgsem R Package. *Structural Equation Modeling: A Multidisciplinary Journal*, 1–28. [DOI](https://doi.org/10.1080/10705511.2026.2617170)

<!--- ## Measurement blocks
The R script "measurement_blocks" describes how to use the Structural-After-Measurement (SAM; Rosseel and Loh, 2024) measurement blocks within the mmgsem package. 
-->
