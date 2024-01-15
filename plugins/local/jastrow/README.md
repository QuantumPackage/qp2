# Jastrow

Information related to the Jastrow factor in trans-correlated calculations.

The main keywords are:
- `j2e_type`
- `j1e_type`
- `env_type`

## j2e_type Options

1. **none:** No 2e-Jastrow is used.

2. **rs-dft:** 2e-Jastrow inspired by Range Separated Density Functional Theory. It has the following shape:
   \begin{equation}
   \tau = \frac{1}{2} \sum_{i,j \neq i} u(\mathbf{r}_i, \mathbf{r}_j),
   \end{equation}
   with,
   \begin{equation}
   u(\mathbf{r}_1, \mathbf{r}_2) = u(r_{12}) = \frac{r_{12}}{2} \left[ 1 - \text{erf}(\mu \, r_{12}) \right] - \frac{\exp\left[- (\mu \, r_{12})^2\right]}{2 \sqrt{\pi} \mu}.
   \end{equation}



## env_type Options

The Jastrow used is multiplied by an envelope \(v\):

\begin{equation}
\tau = \frac{1}{2} \sum_{i,j \neq i} u(\mathbf{r}_i, \mathbf{r}_j) \, v(\mathbf{r}_i) \, v(\mathbf{r}_j)
\end{equation}

- if `env_type` is **none**: No envelope is used.

- if `env_type` is **prod-gauss**: \(v(\mathbf{r}) = \prod_{a} \left(1 - e^{-\alpha_a (\mathbf{r} - \mathbf{R}_a)^2 } \right)\)

- if `env_type` is **sum-gauss**: \(v(\mathbf{r}) = 1 - \sum_{a} \left(1 - c_a e^{-\alpha_a (\mathbf{r} - \mathbf{R}_a)^2 } \right)\)

Here, \(A\) designates the nuclei, and the coefficients and exponents are defined in the tables `enc_coef` and `env_expo` respectively.



## j1e_type Options

The Jastrow used is:

\begin{equation}
\tau = \sum_i u_{1e}(\mathbf{r}_i)
\end{equation}

- if `j1e_type` is **none**: No one-electron Jastrow is used.

- if `j1e_type` is **gauss**: We use \(u_{1e}(\mathbf{r}) = \sum_A \sum_{p_A} c_{p_A} e^{-\alpha_{p_A} (\mathbf{r} - \mathbf{R}_A)^2}\), where the \(c_p\) and \(\alpha_p\) are defined by the tables `j1e_coef` and `j1e_expo`, respectively.

- if `j1e_type` is **charge-harmonizer**: The one-electron Jastrow factor depends on the two-electron Jastrow factor \(u_{2e}\) such that the one-electron term is added to compensate for the unfavorable effect of altering the charge density caused by the two-electron factor:
\begin{equation}
u_{1e}(\mathbf{r}_1) = - \frac{N-1}{2N} \sum_{\sigma} \int d\mathbf{r}_2 \rho^{\sigma}(\mathbf{r}_2) u_{2e}(\mathbf{r}_1, \mathbf{r}_2),
\end{equation}

Feel free to review and let me know if any further adjustments are needed.



