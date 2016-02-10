using Distributions
using GaussianMixtureTest
using KernelEstimator
import PyPlot; const plt = PyPlot

mu_true = [-2.0858,-1.4879]
wi_true = [0.1828,0.8172]
sigmas_true = [0.6735,0.2931]

m = MixtureModel(map((u, v) -> Normal(u, v), mu_true, sigmas_true), wi_true)

T = zeros(100)
P = zeros(100)
for b in 1:100
    srand(b)
    x = rand(m, 1000)
    T[b], P[b] = kstest(x, 2)
    print(b,"->", T[b], " | ")
end

x = rand(m, 1000)
Ttrue = asymptoticdistribution(x, wi_true, mu_true, sigmas_true)

xs = linspace(0.01, 14, 500)
den1=kerneldensity(Ttrue, xeval=xs, kernel=gammakernel, lb=0., h=.2)
#trand = readcsv("/Users/lanfengpan/fa/remote/code/trand2.csv")[:, 1]
den2 = kerneldensity(T, xeval=xs, kernel=gammakernel, lb=0., h=.2)

plt.plot(xs, den1, "k-", xs, den2, "b--")
# write your own tests here
