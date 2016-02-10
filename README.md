# GaussianMixtureTest

[![Build Status](https://travis-ci.org/panlanfeng/GaussianMixtureTest.jl.svg?branch=master)](https://travis-ci.org/panlanfeng/GaussianMixtureTest.jl)

This package implements the EM test for number of components of univariate Gaussian Mixture. The conventional log likelihood test can not be used to test the number of components because the Fisher regularity conditions are violated in Gaussian Mixture case [1].

This package follows the method of [3] but with no regression covariates. Note that the asymptotic distribution of the test statistic is that of the maximum of `C0` dependent $Chi^2(2)$ random variables which has no closed form when the null distribution has more than 1 component. However the p value can be obtained via simulation.

In addition the typical EM algorithm may fail to give a consistent estimate of Gaussian Mixture parameters. This package still uses EM but add a penalty term on the log likelihood which ensures the consistency [2]. 

## Warning

For this method to be valid, a stringent condition (Proposition 6 of [3]) on the true variance is required.

## Usage

To install this package, run

    Pkg.clone("https://github.com/panlanfeng/GaussianMixtureTest.jl.git")

The major functions are `gmm`, `asymptoticdistribution` and `kstest`. `gmm` estimates the parameters via EM algorithm. `asymptoticdistribution` simulates the asymptotic distribution of the test statistic when the number of components is greater than 2. `kstest` conducts the Kasahara-Shimotsu test. 

See also the usage by running

    ?gmm


## Optional Arguments

There are several optional arguments for `gmm`

 - `wi_init` `mu_init` and `sigmas_init`: the initial values
 - `maxiteration`: the number of iterations
 - `tol`: the tolerance of convergence criteria
 - `an` and `sn`: the penalty weight and variance term
 - `taufixed`: for `kstest`, whether fix the `tau` value
 - `whichtosplit` and `tau`: for `kstest`, which component to split and the `split` proportion
 - `mu_lb` and `mu_ub` for `kstest`, the lower and upper limits of components means
    
The optional arguments for `kstest`

 - `vtau`: the finite set of `tau` value
 - `ntrials`: the number of initial values to try
 - `debuginfo`: whether show the debug information

The `npre` of `asymptoticdistribution` means how many random values to generate.

## Examples

See the example code also in [`runtests.jl`](test/runtests.jl)

    using Distributions
    using GaussianMixtureTest

    mu_true = [-2.0858,-1.4879]
    wi_true = [0.0828,0.9172]
    sigmas_true = [0.6735,0.2931]

    m = MixtureModel(map((u, v) -> Normal(u, v), mu_true, sigmas_true), wi_true)
    x = rand(m, 1000);
    
    asymptoticdistribution(x, wi_true, mu_true, sigmas_true, debuginfo=true);
    
    #Estimate the parameters with two components
    wi, mu, sigmas, ml = gmm(x, 2)
    
    #Do the KS test for C=2 v.s. C=3
    kstest(x, 2)

## Acknowledgement

Thanks Dr. Shimotsu and Dr. Kasahara for nicely providing their original R code and their detailed explanations. Several implementation details of this package are borrowed from their R code.

## Reference

 - [1] Chen, J. & Li, P., 2009. Hypothesis Test for Normal Mixture Models: The EM Approach. _the Annals of Statistics_, 37(5 A), pp.2523–2542.

 - [2] Chen, J., Tan, X. & Zhang, R., 2008. Inference for Normal Mixtures in Mean and Variance. _Statistica Sinica_, 18, pp.443–465.
 
 - [3] Kasahara, H. & Shimotsu, K., 2015. Testing the Number of Components in Normal Mixture Regression Models. _Journal of the American Statistical Association_ (to appear), pp.1–33. 
