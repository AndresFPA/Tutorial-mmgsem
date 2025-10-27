# mmgsem Tutorial Script
# Author: Andres Felipe Perez Alonso
# Date: 2025-09-01
# Description: This script describes a multigroup comparison using the mmgsem package
# This script should be read together with the tutorial paper by Perez Alonso, et al. (2025)

# Install the package from github. devtools (or remotes) package is needed.
devtools::install_github("AndresFPA/mmgsem")
# remotes::install_github("AndresFPA/mmgsem") # If devtools is returning errors, you can use remotes to install the package

# Load the library
library("mmgsem")

# Load data
data <- mmgsem::data_example

#### Section 1 - Measurement model ####
# Fit the Step 1 model individually using lavaan
# Define the model using lavaan syntax
S1_model <- '
    # factor loadings
    F1 =~ V1 + V2 + V3 + V4 + V5
    F2 =~ V6 + V7 + V8 + V9 + V10
    F3 =~ V11 + V12 + V13 + V14 + V15
    F4 =~ V16 + V17 + V18 + V19 + V20
'

# Fit the Measurement Model (using a Confirmatory Factor Analysis)
# Note that, at this stage, a measurement invariance testing should be done. We will not go into the detail on this in this tutorial.
# We assume (and know) that metric invariance holds as the data was generated in this way.
MM_fit <- lavaan::cfa(model = S1_model,             # Step 1 model syntax
                      data  = data,                 # Data
                      group = "group",              # Group variables
                      group.equal = "loadings",     # Measurement model constraints (metric invariance)
                      group.partial = c("F1 =~ V2", # MM non-invariances
                                        "F2 =~ V7",
                                        "F3 =~ V12",
                                        "F4 =~ V17"))

#### Section 2 - Do analysis with mmgsem ####
# Define the models of step 1 and step 2 individually
S1_model <- '
    # factor loadings
    F1 =~ V1 + V2 + V3 + V4 + V5
    F2 =~ V6 + V7 + V8 + V9 + V10
    F3 =~ V11 + V12 + V13 + V14 + V15
    F4 =~ V16 + V17 + V18 + V19 + V20
'

S2_model <- '
    # Regression parameters
    F4 ~ F1 + F3
    F3 ~ F1 + F2
'

# Option 1: Use the lavaan object from section 1
mmgsem_fit <- mmgsem(dat     = data,     # Data
                     S1      = S1_model, # Model syntax for step 1 (measurement model)
                     S2      = S2_model, # Model syntax for step 2 (structural model)
                     group   = "group",  # Grouping variable
                     nclus   = 4,        # How many clusters?
                     seed    = 1,        # Seed for replication purposes
                     s1_type = "lavaan", # Type of model used in step 1
                     s1_fit  = MM_fit)   # If we have it, input of step 1 (lavaan object of the measurement model)

summary(mmgsem_fit)

# Option 2: MMGSEM function can run the measurement directly.
# NOTE: It is not possible to use mmgsem() directly to test for measurement invariance. It is recommended to do this separately using lavaan
# mmgsem_fit <- mmgsem::mmgsem(dat    = data,      # Data
#                              S1     = S1_model,  # Model syntax for step 1 (measurement model)
#                              S2     = S2_model,  # Model syntax for step 2 (structural model)
#                              group  = "group",   # Grouping variable
#                              nclus  = 4,         # How many clusters?
#                              seed   = 1,         # Seed for replication purposes
#                              s1_type = "lavaan", # Type of model used in step 1
#                              s1_fit = NULL)      # Automatically set to NULL if we do not have a lavaan object for the measurement model
#
# mmgsem::summary(mmgsem_fit)

#### Section 3 - Model Selection ####
# In empirical research, we do not know the true number of clusters for the data at hand. Thus, a 'model selection' procedure
# is often done to find the most appropriate number of clusters. For this, one runs several models (e.g., from 1 to 6 clusters) and
# selects the best one based on an specific criteria (e.g., BIC)
# In mmgsem, one can use the ModelSelection() function, which has a similar set of arguments.
# Fit several models simultaneously:
set.seed(1) #set seed to replicate results
modelSelection_fit <- ModelSelection(dat      = data,     # Data
                                     S1       = S1_model, # Model syntax for step 1 (measurement model)
                                     S2       = S2_model, # Model syntax for step 2 (structural model)
                                     group    = "group",  # Grouping variable
                                     nclus    = c(1,6),   # Instead of a specific number of clusters, set the lower and upper limit of the models we want to run (e.g., from 1 to 6 clusters)
                                     s1_type  = "lavaan", # Type of model used in step 1, 
                                     seed     = 1,        # Seed for replication purposes
                                     s1_fit   = MM_fit)   # If we have it, input of step 1 (lavaan object of the measurement model)

summary(model = modelSelection_fit, model_selection = T)

# Check the results from some model selection critera in scree plots
# Convex Hull and AIC
CH_plot   <- plot(ModelSel = modelSelection_fit, criteria = "Chull")
AIC3_plot <- plot(ModelSel = modelSelection_fit, criteria = "AIC3")
cowplot::plot_grid(CH_plot$Chull, AIC3_plot$AIC3, labels = c("Convex Hull", "AIC_3"),
                   hjust = c(-2.8, -6), ncol = 1)

# Given that most model selection measures select the 4-cluster model, we decide to keep that one
# One can easily extract the selected model
selected_fit <- extract(object = modelSelection_fit, nclus = 4)
summary(selected_fit)

# Check posterior matrix to see the cluster memberships
print(selected_fit$posteriors)

#### Section 4 - Standard Errors and Hypothesis Testing ####
# As the standard error computation can take a significant amount of time. It is not integrated in the mmgsem() function

# To compute the standard errors, we use the compute_se() function
# standard_errors <- compute_se(object = selected_fit) # Commented out due to computation time

# Alternatively, load the pre-computed standard errors for this tutorial (or use naive = T to reduce the computation time)
standard_errors <- mmgsem::standard_errors

## Hypothesis testing 1: Are the betas different from 0?
# To check the se in a clearer way, one can use the standard_errors object in the summary function
summary(model = selected_fit, se = standard_errors)

## Hypothesis testing 2: Are the betas different across clusters?
# Wald test to check if there is at least one parameter that differs across clusters
test.mmgsem(model = selected_fit, se = standard_errors)

# For a more detailed testing, we can check differences in parameters across all pairs of clusters
test.mmgsem(model = selected_fit, se = standard_errors, multiple_comparison = TRUE)

#### Section 4 - Using hypothesis testing to support model selection ####
# As an example, check the 5-cluster model
five_K_fit <- extract(object = modelSelection_fit, nclus = 5)
summary(model = five_K_fit)

# Compute the corrected standard errors (takes a long time)
# five_K_se <- compute_se(object = five_K_fit) # Commented out due to computation time

# Alternatively, load the pre-computed standard errors for this tutorial (or use naive = T to reduce the computation time)
five_k_se <- mmgsem::five_K_se

# Check hypothesis testing
summary(model = five_K_fit, se = five_K_se) # Group 21 is isolated in cluster 4 (and should actually be part of cluster 5)

# Multiple comparison testing
test.mmgsem(model = five_K_fit, se = five_K_se, multiple_comparison = TRUE) # Clusters 4 and 5 never significantly differ 




