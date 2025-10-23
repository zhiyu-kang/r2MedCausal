# r2MedCausal

`r2MedCausal` provides methods for estimating high-dimensional total mediation effect in case-control studies using an R-squared based measure. This approach accounts for ascertainment bias, weak mediators, and bidirectional effects using a cross-fitted, modified Haseman-Elston regression-based estimation procedure.

---

## Installation

You can install `r2MedCausal` package using `devtools`:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install r2MedCausal from GitHub
devtools::install_github("zhiyu-kang/r2MedCausal")
library(r2MedCausal)
```

---

## Usage Example

Here’s a simple example to generate synthetic data and run mediation analysis:

```r
# Load the package
library(r2MedCausal)

# Generate synthetic data
data <- data_generate(
  p = 2000, N = 500, K = 0.05, r = 0.5, res.var = 0.2,
  pi_11 = 0.09, pi_alpha = 0.3, pi_beta = 0.3,
  sig2_a = 0.3, sig2_11 = 0.5, sig2_01 = 0.3, seed = 23
)

# Extract variables
X <- data$x
M <- data$M
Y <- data$y
K <- 0.05

# Run mediation analysis
result <- r2_estimation_cf(X, M, Y, K)
print(result)
```

You can compute confidence intervals using jackknife resampling. 

```r
# Set the number of cores for parallel computation
num_cores <- 6  # Adjust based on your system capabilities

# Compute confidence intervals
ci_result <- confidence_interval(num_cores, X, M, Y, K)
print(ci_result)
```

 **Warning:**  
- Jackknife resampling involves removing one observation at a time and recomputing the mediation effect.  
- This makes the computation **significantly slower** than a single estimation.  
- For large datasets, consider using a high-performance computing cluster.


## Citing r2MedCausal

If you use `r2MedCausal` in your research, please cite the following paper:

Zhiyu Kang, Li Chen, Peng Wei, Zhichao Xu, Chunlin Li, Tianzhong Yang (2025). **Estimation of total mediation effect for a binary trait in a case-control study for high-dimensional omics mediators.** *bioRxiv*. https://doi.org/10.1101/2025.01.28.635396
