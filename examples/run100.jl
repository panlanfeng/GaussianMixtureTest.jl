using GaussianMixtureTest
using Distributions
using KernelEstimator
using RCall

mu_true = [-2.0858,-1.4879]
wi_true = [0.0828,0.9172]
sigmas_true = [0.6735,0.2931]

m = MixtureModel(map((u, v) -> Normal(u, v), mu_true, sigmas_true), wi_true)

T1 = zeros(100)
P = zeros(100)
for b in 1:100
    srand(b)
    x = rand(m, 500)
    T1[b], P[b] = kstest(x, 2)
    print(b,"->", T1[b], " | ")
end

x = rand(m, 500)
Ttrue = asymptoticdistribution(x, wi_true, mu_true, sigmas_true)

xs = linspace(0.01, 14, 500)
den1=kerneldensity(Ttrue, xeval=xs, kernel=gammakernel, lb=0.)
den2 = kerneldensity(T1, xeval=xs, kernel=gammakernel, lb=0.)


@rput xs den1 den2 T1 P
rprint(""" 
hist(T1, breaks=15, freq=F, ylim=c(0, .4))
lines(xs, den1, lwd=2)
lines(xs, den2, lwd=2, col="blue")
""")

rprint("""
z = qnorm(1-P)
hist(z, freq=F, xlim=c(-4,4))
lines(density(z))
curve(dnorm, -4, 4, col="red",add=T)
NULL
""")


## The fitting

import GaussianMixtureTest
using Distributions
using KernelEstimator
using RCall

mu_true = [-2.0858,-1.4879]
wi_true = [0.0828,0.9172]
sigmas_true = [0.6735,0.2931]

m = MixtureModel(map((u, v) -> Normal(u, v), mu_true, sigmas_true), wi_true)

x = rand(m, 500)
wi, mu, sigmas = GaussianMixtureTest.gmm(x, 2)
mhat = MixtureModel(map((u, v) -> Normal(u, v), mu, sigmas), wi)

xs = linspace(-5, 2, 500)
dentrue = pdf(m, xs)
denhat = pdf(mhat, xs)

@rput xs dentrue denhat
rprint(""" 
plot(xs, dentrue, lwd=2, type="l")
lines(xs, denhat, lwd=2, col="blue")
""")

Ttrue = GaussianMixtureTest.asymptoticdistribution(x, wi_true, mu_true, sigmas_true)
xs = linspace(0.01, 14, 500)
den1=kerneldensity(Ttrue, xeval=xs, kernel=gammakernel, lb=0.)

@rput xs den1
rprint(""" 
plot(xs, den1, type="l")
""")
