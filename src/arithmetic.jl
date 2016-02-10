
function log1pexp!(res::AbstractArray{Float64}, x::AbstractArray{Float64}, res2::AbstractArray{Float64}, n::Int64=length(x))
    
    copy!(res2, x)
    Yeppp.exp!(res, x)
    log1p!(res, res, n)
    for i in 1:n
        @inbounds res[i]=mylog1pexp(res[i], res2[i])
    end
    # for i in 1:n
    #     @inbounds res[i]=log1pexp(x[i])
    # end
end
mylog1pexp(x::AbstractFloat, x2::AbstractFloat)=x2 < 18.0 ? x : x2 < 33.3 ? x2+exp(-x2) : x2
function log1pexp!(res::AbstractArray{Float64}, x::AbstractArray{Float64}, n::Int64=length(x))
     
    res2=copy(x)
    log1pexp!(res, x, res2, n)
end
# log1pexp(x::AbstractFloat) = x < 18.0 ? log1p(exp(x)) : x < 33.3 ? x + exp(-x) : x
log1pexp!(x::AbstractArray{Float64}, n::Int64=length(x))=log1pexp!(x, x, n)

function sumexp{T<:Real}(x::AbstractArray{T})
    isempty(x) && return -Inf
    u = maximum(x)
    s = 0.
    for i in 1:length(x)
        @inbounds s += exp(x[i] - u)
    end
    s * exp(u)
end
function sumexp{T<:Real}(x::AbstractArray{T}, coef::AbstractArray{T})
    isempty(x) && return -Inf
    u = maximum(x)
    s = 0.
    for i in 1:length(x)
        @inbounds s += coef[i]*exp(x[i] - u)
    end
    s * exp(u)
end

function ratiosumexp!{T<:Real}(x::AbstractArray{T}, coef::AbstractArray{T}, s::AbstractArray{T}, ncomponent::Int)
    #length(x) != length(coef) && error("Length should be the same!")
    #isempty(x) && return -Inf
    u = maximum(x)
    for i in 1:ncomponent
        @inbounds s[i] = coef[i]*exp(x[i] - u)
    end
    divide!(s, s, sum(s), ncomponent)
    #s ./ sum(s) #* exp(u)
    nothing
end
function ratiosumexp!{T<:Real}(x::AbstractArray{T}, coef::AbstractArray{T}, s::AbstractMatrix{T}, irow::Int, ncomponent::Int)

    u = maximum(x)
    ssum=0.0
    @inbounds for i in 1:ncomponent
        tmp=coef[i]*exp(x[i] - u)
        s[irow, i] = tmp
        ssum += tmp
    end
    #ssum = sum(s[irow, :])
    for i in 1:ncomponent
        @inbounds s[irow, i] = s[irow, i] / ssum
    end
    nothing
end
function add!(res::AbstractArray{Float64}, x::AbstractArray{Float64}, y::Float64, n::Int64=length(x))
   for i in 1:n
       @inbounds res[i] = x[i] + y
   end
   nothing
end

add!(res::AbstractArray{Float64}, x::Float64, y::AbstractArray{Float64}, n::Int64=length(y)) = add!(res, y, x, n)

add!(x::AbstractArray{Float64}, y::Float64, n::Int64=length(x))=add!(x, x, y, n)
add!(x::Float64, y::AbstractArray{Float64}, n::Int64=length(y))=add!(y, y, x, n)

plusone!(res::AbstractArray{Float64}, x::AbstractArray{Float64}, n::Int64=length(x)) = add!(res, x, 1.0, n)
plusone!(x::AbstractArray{Float64}, n=length(x)) = plusone!(x, x, n)


function divide!(res::AbstractArray{Float64}, x::AbstractArray{Float64}, y::AbstractArray{Float64}, n::Int64=length(x))
    for i in 1:n
        @inbounds res[i] = x[i] / y[i]
    end
    nothing
end

divide!(res::AbstractArray{Float64}, x::AbstractArray{Float64}, y::Float64, n::Int64=length(x)) = multiply!(res, x, 1/y, n)
function divide!(res::AbstractArray{Float64}, x::Float64, y::AbstractArray{Float64}, n::Int64=length(y))
    for i in 1:n
        @inbounds res[i] = x / y[i]
    end
    nothing
end
divide!(x::AbstractArray{Float64}, y::AbstractArray{Float64}, n::Int64=length(x)) = divide!(x, x, y, n)
# divide!(x::Float64, y::AbstractArray{Float64}, n::Int64=length(y)) = divide!(y, x, y, n)
divide!(x::AbstractArray{Float64}, y::Float64, n::Int64=length(x)) = divide!(x, x, y, n)



rcp!(res::AbstractArray{Float64}, x::AbstractArray{Float64}, n::Int64=length(x)) = divide!(res, 1.0, x, n)
rcp!(x::AbstractArray{Float64}, n::Int64=length(x))=rcp!(x, x, n)



function multiply!(res::AbstractArray{Float64}, x::AbstractArray{Float64}, y::Float64, n::Int64=length(x))
    for i in 1:n
        @inbounds res[i] = x[i] * y
    end
    nothing
end

function multiply!(res::AbstractArray{Float64}, x::Float64, y::AbstractArray{Float64}, n::Int64=length(y))
    for i in 1:n
        @inbounds res[i] = x * y[i]
    end
    nothing
end

# multiply!(x::Float64, y::AbstractArray{Float64}, n::Int64=length(y)) = multiply!(y, x, y, n)
multiply!(x::AbstractArray{Float64}, y::Float64, n::Int64=length(x)) = multiply!(x, x, y, n)

function negate!(res::AbstractArray{Float64}, x::AbstractArray{Float64}, n::Int64=length(x))
   for i in 1:n
       @inbounds res[i] = -x[i]
   end
   nothing
end
negate!(x::AbstractArray{Float64}, n::Int64=length(x)) = negate!(x, x, n)

function negateiffalse!(x::AbstractArray{Float64}, y::AbstractArray{Bool, 1}, n::Int64=length(x))
    for i in 1:n
        @inbounds x[i] = ifelse(y[i], x[i], -x[i])
    end
    nothing
end

function negateiftrue!(x::AbstractArray{Float64}, y::AbstractArray{Bool, 1}, n::Int64=length(x))
    for i in 1:n
        @inbounds x[i] = ifelse(y[i], -x[i], x[i])
    end
    nothing
end

function log1p!(res::AbstractArray{Float64}, x::AbstractArray{Float64}, n::Int64=length(x))
    add!(res, x, 1., n)
    Yeppp.log!(res, res)
    res
end
log1p!(x::AbstractArray{Float64}, n::Int64=length(x)) = log1p!(x, x, n)

# -log(1+exp(-xy))
function loglogistic!(x::AbstractArray{Float64}, y::AbstractArray{Bool, 1}, n::Int=length(x))
    negateiftrue!(x, y, n)
    log1pexp!(x, x, n)
    negate!(x, x, n)
    nothing
end

# 1/(1+exp(-x))
function logistic!(res::AbstractArray{Float64}, x::AbstractArray{Float64}, n::Int=length(x))
    
    negate!(res, x, n)
    Yeppp.exp!(res, res)
    add!(res, res, 1.0, n)
    rcp!(res, res, n)
    res
end
logistic!(x::AbstractArray{Float64}, n::Int=length(x)) = logistic!(x, x, n)


function relocate!(res::AbstractArray{Float64}, ga::Vector{Float64}, groupindex::IntegerVector, N::Int)
    for i in 1:N
        @inbounds res[i] = ga[groupindex[i]]
    end
    nothing
end

function sumby!(r::AbstractArray, y::AbstractArray, x::IntegerArray, levels::IntUnitRange)
	k = length(levels)
	length(r) == k || raise_dimerror()

	m0 = levels[1]
	m1 = levels[end]
	b = m0 - 1

	@inbounds for i in 1 : length(x)
		xi = x[i]
		if m0 <= xi <= m1
			r[xi - b] += y[i]
		end
	end
	#return r
end
function sumsqby!(r::AbstractArray, y::AbstractArray, x::IntegerArray, levels::IntUnitRange)
	k = length(levels)
	length(r) == k || raise_dimerror()

	m0 = levels[1]
	m1 = levels[end]
	b = m0 - 1

	@inbounds for i in 1 : length(x)
		xi = x[i]
		if m0 <= xi <= m1
			r[xi - b] += abs2(y[i])
		end
	end
	#return r
end

function H1(y, mu, sigmas)
    (y .- mu)./sigmas./sigmas
end
function H2(y, mu, sigmas)
    z = (y .- mu)./sigmas
    (z.^2 .-1)./sigmas^2./2
end

function H3(y, mu, sigmas)
    z = (y .- mu)./sigmas
    (z.^3 .-3.*z) ./ sigmas^3./6
end

function H4(y, mu, sigmas)
    z = (y .- mu)./sigmas
    (z.^4 .-6.*z.^2 .+ 3) ./ sigmas^4 ./ 24
end
