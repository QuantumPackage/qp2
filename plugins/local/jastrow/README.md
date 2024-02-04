# Jastrow

Information related to the Jastrow factor in trans-correlated calculations.

The main keywords are:
- `j2e_type`
- `j1e_type`
- `env_type`

## j2e_type Options

1. **None:** No 2e-Jastrow is used.

2. **Mu:** 2e-Jastrow inspired by Range Separated Density Functional Theory. It has the following shape:
   <p align="center">
      <img src="https://latex.codecogs.com/png.image?%5Cinline%20%5Clarge%20%5Cdpi%7B200%7D%5Cbg%7Bwhite%7D%5Ctau=%5Cfrac%7B1%7D%7B2%7D%5Csum_%7Bi,j%5Cneq%20i%7Du(%5Cmathbf%7Br%7D_i,%5Cmathbf%7Br%7D_j)">
   </p>
   with,
   <p align="center">
   <img src="https://latex.codecogs.com/png.image?%5Cinline%20%5Clarge%20%5Cdpi%7B200%7D%5Cbg%7Bwhite%7D%20u(%5Cmathbf%7Br%7D_1,%5Cmathbf%7Br%7D_2)=u(r_%7B12%7D)=%5Cfrac%7Br_%7B12%7D%7D%7B2%7D%5Cleft%5B1-%5Ctext%7Berf%7D(%5Cmu%20r_%7B12%7D)%5Cright%5D-%5Cfrac%7B%5Cexp%5B-(%5Cmu%20r_%7B12%7D)%5E2%5D%7D%7B2%5Csqrt%7B%5Cpi%7D%5Cmu%7D">
   </p>

3. **Mu_Nu:** A valence and a core correlation terms are used
   <p align="center">
      <img src="https://latex.codecogs.com/png.image?\inline&space;\huge&space;\dpi{110}\bg{white}&space;u(\mathbf{r}_1,\mathbf{r}_2)=u(\mu;r_{12})\,v(\mathbf{r}_1)\,v(\mathbf{r}_2)&plus;u(\nu;r_{12})[1-v(\mathbf{r}_1)\,v(\mathbf{r}_2)]">
   </p>
   with envelop \(v\).


## env_type Options

The 2-electron Jastrow is multiplied by an envelope \(v\):
<p align="center">
      <img src="https://latex.codecogs.com/png.image?%5Cinline%20%5Clarge%20%5Cdpi%7B200%7D%5Cbg%7Bwhite%7D%5Ctau=%5Cfrac%7B1%7D%7B2%7D%5Csum_%7Bi,j%5Cneq%20i%7Du(%5Cmathbf%7Br%7D_i,%5Cmathbf%7Br%7D_j)%5C,v(%5Cmathbf%7Br%7D_i)%5C,v(%5Cmathbf%7Br%7D_j)">
</p>

- if `env_type` is **None**: No envelope is used.

- if `env_type` is **Prod_Gauss**:
  <p align="center">
     <img src="https://latex.codecogs.com/png.image?%5Cinline%20%5Clarge%20%5Cdpi%7B200%7D%5Cbg%7Bwhite%7D%20v(%5Cmathbf%7Br%7D)=%5Cprod_%7BA%7D%5Cleft(1-e%5E%7B-%5Calpha_A(%5Cmathbf%7Br%7D-%5Cmathbf%7BR%7D_A)%5E2%7D%5Cright)">
   </p>

- if `env_type` is **Sum_Gauss**:
  <p align="center">
     <img src="https://latex.codecogs.com/png.image?%5Cinline%20%5Clarge%20%5Cdpi%7B200%7D%5Cbg%7Bwhite%7D%20v(%5Cmathbf%7Br%7D)=1-%5Csum_%7BA%7Dc_A%20e%5E%7B-%5Calpha_A(%5Cmathbf%7Br%7D-%5Cmathbf%7BR%7D_A)%5E2%7D">
  </p>

Here, \(A\) designates the nuclei, and the coefficients and exponents are defined in the tables `env_coef` and `env_expo` respectively.


## j1e_type Options

The 1-electron Jastrow used is:
<p align="center">
   <img src="https://latex.codecogs.com/png.image?%5Cinline%20%5Clarge%20%5Cdpi%7B200%7D%5Cbg%7Bwhite%7D%5Ctau=%5Csum_i%20u_%7B1e%7D(%5Cmathbf%7Br%7D_i)">
</p>

- if `j1e_type` is **None**: No one-electron Jastrow is used.

- if `j1e_type` is **Gauss**: We use
<p align="center">
   <img src="https://latex.codecogs.com/png.image?%5Cinline%20%5Clarge%20%5Cdpi%7B200%7D%5Cbg%7Bwhite%7Du_%7B1e%7D(%5Cmathbf%7Br%7D)=%5Csum_A%5Csum_%7Bp_A%7Dc_%7Bp_A%7De%5E%7B-%5Calpha_%7Bp_A%7D(%5Cmathbf%7Br%7D-%5Cmathbf%7BR%7D_A)%5E2%7D">
</p>
<img src="https://latex.codecogs.com/png.image?%5Cinline%20%5Clarge%20%5Cdpi%7B200%7D%5Cbg%7Bwhite%7D%20c_%7Bp_A%7D%5C,%5Ctext%7Band%7D%5C,%5Calpha_%7Bp_A%7D"> 

are defined by the tables `j1e_coef` and `j1e_expo`, respectively.

- if `j1e_type` is **Charge_Harmonizer**: The one-electron Jastrow factor aims to offset the adverse impact of modifying the charge density induced by the two-electron factor
  <p align="center">
     <img src="https://latex.codecogs.com/png.image?%5Cinline%20%5Clarge%20%5Cdpi%7B200%7D%5Cbg%7Bwhite%7Du_%7B1e%7D(%5Cmathbf%7Br%7D_1)=-%5Cfrac%7BN-1%7D%7B2N%7D%5C,%5Csum_%7B%5Csigma%7D%5C,%5Cint%20d%5Cmathbf%7Br%7D_2%5C,%5Crho%5E%7B%5Csigma%7D(%5Cmathbf%7Br%7D_2)%5C,u_%7B2e%7D(%5Cmathbf%7Br%7D_1,%5Cmathbf%7Br%7D_2)">
  </p>

- if `j1e_type` is **Charge_Harmonizer_AO**: The one-electron Jastrow factor **Charge_Harmonizer** is fitted by the product of atomic orbitals:
  <p align="center">
     <img src="https://latex.codecogs.com/png.image?\inline&space;\large&space;\dpi{300}\bg{white}&space;u_{1e}(\mathbf{r})=\sum_{\alpha,\beta}C_{\alpha,\beta}\chi_{\alpha}(\mathbf{r})\chi_{\beta}(\mathbf{r})">
  </p>

