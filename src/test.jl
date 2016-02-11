"""
    asymptoticdistribution(x, wi, mu, sigmas);

Simulate the asymptotic distribution of test statistic for `kstest`. 
`nrep` is the number of random values to generate.
`debuginfo` is whether to show the debug information.
When number of components of null distribution is greater than 1, the test statistic has no closed form asymptotic distribution. When the null distribution is just normal, the asymptotic distribution is just `Chisq(2)`.

"""
function asymptoticdistribution(x::RealVector{Float64}, wi::Vector{Float64}, mu::Vector{Float64}, sigmas::Vector{Float64}; nrep::Int=10000, debuginfo::Bool=false)

    n = length(x)
    C = length(wi)
    m = MixtureModel(map((u, v) -> Normal(u, v), mu, sigmas), wi)
    llC = zeros(n, C)
    S_π = zeros(n, C-1)
    S_μσ = zeros(n, 2*C)
    S_λ = zeros(n, 2*C)
    ll = logpdf(m, x)
    for i in 1:n, kcom in 1:C
        llC[i, kcom] = logpdf(m.components[kcom], x[i])
    end

    for kcom in 1:(C-1)
        S_π[:, kcom] = exp(llC[:, kcom] .- ll) .- exp(llC[:, C] .- ll) #(llC[:, kcom] .- llC[:, C]) ./ ll
    end
    for i in 1:n
        for kcom in 1:C
            llC[i, kcom] = exp(log(wi[kcom]) + llC[i, kcom] - ll[i])
            S_μσ[i, 2*kcom-1] = H1(x[i], mu[kcom], sigmas[kcom]) * llC[i, kcom]
            S_μσ[i, 2*kcom] = H2(x[i], mu[kcom], sigmas[kcom]) * llC[i, kcom]
            S_λ[i, 2*kcom-1] = H3(x[i], mu[kcom], sigmas[kcom]) * llC[i, kcom]
            S_λ[i, 2*kcom] = H4(x[i], mu[kcom], sigmas[kcom]) * llC[i, kcom]
        end
    end
    S_η = hcat(S_π, S_μσ)
    debuginfo && println(round(llC[1:5,:], 6))
    debuginfo && println(sum(S_η, 1))
    debuginfo && println(sum(S_λ, 1))
    I_η = S_η'*S_η./n
    I_λη = S_λ'*S_η./n
    I_λ = S_λ'*S_λ./n
    I_all = vcat(hcat(I_η, I_λη'), hcat(I_λη, I_λ))
    if 1/cond(I_all) < eps(Float64)
        D, V = eig(I_all)
        debuginfo && println(D)
        tol2 = maximum(abs(D)) * 1e-14
        D[D.<tol2] = tol2
        I_all = V*diagm(D)*V'
    end
    debuginfo && println(round(cor(S_η), 6))
    debuginfo && println(round(cor(S_λ), 6))
    I_λ_η = I_all[(3*C):(5*C-1), (3*C):(5*C-1)] - I_all[(3*C):(5*C-1), 1:(3*C-1)] * inv(I_all[1:(3*C-1), 1:(3*C-1)]) * I_all[1:(3*C-1),(3*C):(5*C-1)]
    debuginfo && println(round(I_λ_η, 6))
    #I_λ_η=(I_λ_η .+ I_λ_η')./2
    D, V = eig(I_λ_η)
    D[D.<0.] = 0.
    debuginfo && println(D)
    I_λ_η2 = V * diagm(sqrt(D)) * V'
    u = randn(nrep, 2*C) * I_λ_η2
    EM = zeros(nrep, C)
    T = zeros(nrep)
    for kcom in 1:C
        EM[:, kcom] = sum(u[:, (2*kcom-1):(2*kcom)] * inv(I_λ_η[(2*kcom-1):(2*kcom), (2*kcom-1):(2*kcom)]) .* u[:, (2*kcom-1):(2*kcom)], 2)
    end
    for i in 1:nrep
        T[i] = maximum(EM[i, :])
    end
    debuginfo && println(EM[1:10,:])
    T
end

function maxll_ctau(x::RealVector, C1::Int, wi_C1::Vector{Float64}, whichtosplit::Int,
   tau::Real, mu_lb::Vector{Float64}, mu_ub::Vector{Float64}, sigmas_lb::Vector{Float64}, sigmas_ub::Vector{Float64},
   ntrials::Int, sn::Vector{Float64}, an::Real, debuginfo::Bool, tol::Real)
   
    n = length(x)
    tau = min(tau, 1-tau)

    wi = repmat(wi_C1, 1, 4*ntrials)
    mu = zeros(C1, 4*ntrials)
    sigmas = ones(C1, 4*ntrials)
    ml = -Inf .* ones(4*ntrials)
    for i in 1:4*ntrials
        mu[:, i] = rand(C1) .* (mu_ub .- mu_lb) .+ mu_lb
        sigmas[:, i] = rand(C1) .* (sigmas_ub .- sigmas_lb) .+ sigmas_lb

        wi[:, i], mu[:, i], sigmas[:, i], ml[i] =
             gmm(x, C1, wi[:, i], mu[:, i], sigmas[:, i],
             whichtosplit=whichtosplit, tau=tau,
             mu_lb=mu_lb, mu_ub=mu_ub,
             maxiteration=100, sn=sn, an=an,
             taufixed=true, tol=tol)
    end

    mlperm = sortperm(ml)
    for j in 1:ntrials
        i = mlperm[4*ntrials+1 - j] # start from largest ml
        wi[:, i], mu[:, i], sigmas[:, i], ml[i] =
            gmm(x, C1, wi[:, i], mu[:, i], sigmas[:, i],
            whichtosplit=whichtosplit, tau=tau,
            mu_lb=mu_lb, mu_ub=mu_ub,
            sn=sn, an=an, taufixed=true,
            tol=tol)
    end

    mlmax, imax = findmax(ml[mlperm[(3*ntrials+1):4*ntrials]])
    imax = mlperm[3*ntrials+imax]

    re=gmm(x, C1, wi[:, imax], mu[:, imax], sigmas[:, imax],
         maxiteration=3, an=an, sn=sn, tol=0.)
    debuginfo && println("Trial:", re)
    return(re[4])
    
end
 
"""
    kstest(x, C0)
    
Do the EM test under null Hypothesis of `C0` components.
If rejected, then it suggest the true number of components is greater than `C0`.
Optional arguments for `kstest`

 - `vtau`: the finite set of `tau` value
 - `ntrials`: the number of initial values to try
 - `debuginfo`: whether show the debug information

"""
function kstest(x::RealVector{Float64}, C0::Int; vtau::Vector{Float64}=[.5,.3,.1;],
    ntrials::Int=25, debuginfo::Bool=false, tol::Real=0.001)
    
    C1 = C0+1
    n = length(x)

    wi_init, mu_init, sigmas_init, ml_C0 = gmm(x, C0)

    debuginfo && println(wi_init, mu_init, sigmas_init, ml_C0)
    if C0 > 1
        trand=GaussianMixtureTest.asymptoticdistribution(x, wi_init, mu_init, sigmas_init)
    end

    if debuginfo
        println("ml_C0=", ml_C0)
    end
    minx = minimum(x)
    maxx = maximum(x)

    or = sortperm(mu_init)
    wi0 = wi_init[or]
    mu0 = mu_init[or]
    sigmas0 = sigmas_init[or]

    an = decidepenalty(wi0, mu0, sigmas0, n)
    lr = 0.0
    lrv = zeros(length(vtau), C0)
    for whichtosplit in 1:C0, i in 1:length(vtau)
        ind = [1:whichtosplit, whichtosplit:C0;]
        if C1==2
            mu_lb = minx .* ones(2)
            mu_ub = maxx .* ones(2)
        elseif C1>2
            mu_lb = [minx, (mu0[1:(C0-1)] .+ mu0[2:C0])./2;]
            mu_ub = [(mu0[1:(C0-1)] .+ mu0[2:C0])./2, maxx;]
            mu_lb = mu_lb[ind]
            mu_ub = mu_ub[ind]
        end
        sigmas_lb = 0.25 .* sigmas0[ind]
        sigmas_ub = 2 .* sigmas0[ind]

        wi_C1 = wi0[ind]
        wi_C1[whichtosplit] = wi_C1[whichtosplit]*vtau[i]
        wi_C1[whichtosplit+1] = wi_C1[whichtosplit+1]*(1-vtau[i])

        lrv[i, whichtosplit] = maxll_ctau(x, C1, wi_C1, whichtosplit,
           vtau[i], mu_lb, mu_ub, sigmas_lb, sigmas_ub,
           ntrials, sigmas0[ind], an, debuginfo, tol)
       if debuginfo
           println(whichtosplit, " ", vtau[i], "->",
           lrv[i, whichtosplit])
       end
   end
   lr = maximum(lrv)
   if debuginfo
       println("lr=", lr)
   end
   Tvalue = 2*(lr - ml_C0)
   if C0 == 1
       pvalue = 1 - cdf(Chisq(2), Tvalue)
   else
       pvalue = mean(trand .> Tvalue)
   end
   return(Tvalue, pvalue)
end
