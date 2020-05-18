
module RationalPolynomials
    using Polynomials
    include("PolynomialUtils.jl")
    

    export RationalPolynomial, fromrootsandpoles, rationalFit, reduce
    
    struct RationalPolynomial{T <: Number}
        numerator::Poly{<:T}
        denominator::Poly{<:T}
        function RationalPolynomial(numcoeffs::AbstractArray{<:T}, dencoeffs::AbstractArray{<:T}, var::Symbol = :x) where {T <: Number}
            return new{T}(Poly(numcoeffs, var), Poly(dencoeffs, var))
        end

        function RationalPolynomial(numerator::Poly{T}, denominator::Poly{T}) where {T <: Number}
            return new{T}(numerator, denominator)
        end
    end

    (p::RationalPolynomial{T})(x::T) where {T} =  p.numerator(x) / p.denominator(x)

    (p::RationalPolynomial{T})(x::AbstractArray{<:T}) where {T} = p.numerator(x) ./ p.denominator(x)

    function fromrootsandpoles(r::AbstractArray{T}, p::AbstractArray{T}; var::Symbol = :x) where {T <: Number}
        return RationalPolynomial(poly(r, var), poly(p, var))
    end

    function setRationalEquation(degN::Unsigned, degD::Unsigned, x::Number, y::Number)
        nums = [x^i for i in 0:degN]
        dens = [x^i for i in 1:degD]
        dens *= -y
        [nums; dens], y
    end

    function setRationalEquationMatrix(degN::Unsigned, degD::Unsigned, xys::AbstractArray{<:Tuple{<:Number, <:Number}})
        l = length(xys)
        v = Array{Number}(undef, 0)
        as = Array{Number}(undef, 0)
        for (x, y) in xys
            bs, c = setRationalEquation(degN, degD, x, y)
            as = vcat(as, bs)
            v = vcat(v, c)
        end
        transpose(reshape(as, (l, l))), v
    end

    function decomposeToRational(v::AbstractArray{<:Number}, degN::Unsigned, degD::Unsigned)
        vn = v[1:(degN  + 1)]
        vd = vcat([1], v[(degN + 2):end])
        RationalPolynomial(vn, vd)
    end

    function rationalFit(xys::AbstractArray{Tuple{<:Number,<:Number}}, degN::Unsigned, degD::Unsigned)
        if length(xys) != degN + degD + 1
            throw(DimensionMismatch("Sum of degrees mismatches the number of given points"))
        end
        A, v = setRationalEquationMatrix(degN, degD, xys)
        decomposeToRational(A\v, degN, degD)
    end

    function reduce(r::RationalPolynomial{T}) where {T}
        gcdp = gcd(r.numerator, r.denominator)
        return RationalPolynomial(PolyUtils.polyLongDivision(r.numerator, gcdp), PolyUtils.polyLongDivision(r.denominator, gcdp))
    end

    inverse(r::RationalPolynomial{T}) where {T} =  RationalPolynomial(r.denominator, r.numerator)

    function Base.:+(r1::RationalPolynomial{T}, r2::RationalPolynomial{S}) where {T, S}
        return reduce(RationalPolynomial(
            r1.numerator * r2.denominator + r1.denominator * r2.numerator, 
            r1.denominator * r2.denominator))
    end

    Base.:+(r::RationalPolynomial{T}, c::S) where {T, S<:Number} = RationalPolynomial(r.numerator + c * r.denominator, r.denominator)

    Base.:+(c::S, r::RationalPolynomial{T}) where {T, S<:Number} = r + c

    Base.:+(r::RationalPolynomial{T}, p::Poly{S}) where {T, S} = reduce(RationalPolynomial(r.numerator + p * r.denominator, r.denominator))

    Base.:+(p::Poly{S}, r::RationalPolynomial{T}) where {T, S} = r + p
    
    function Base.:*(r1::RationalPolynomial{T}, r2::RationalPolynomial{S}) where {T,S}
        return RationalPolynomial(r1.numerator * r2.numerator, r1.denominator * r2.denominator)
    end

    Base.:*(r1::RationalPolynomial{T}, c::S) where {T, S<:Number} = RationalPolynomial(r1.numerator * c, r1.denominator)

    Base.:*(c::S, r::RationalPolynomial{T}) where {T, S<:Number} = r * c

    Base.:*(r::RationalPolynomial{T}, p::Poly{S}) where {T, S} = reduce(RationalPolynomial(r.numerator * p, r.denominator))

    Base.:*(p::Poly{S}, r::RationalPolynomial{T}) where {T, S} = r * p

    Base.:-(c::S, r::RationalPolynomial{T}) where {T, S<:Number} = c + r * -1

    Base.:-(r::RationalPolynomial{T}, c::S) where {T, S<:Number} = r + -1 * c

    Base.:-(r1::RationalPolynomial{T}, r2::RationalPolynomial{S}) where {T, S} = r1 + (-1 * r2)

    Base.:-(r::RationalPolynomial{T}, p::Poly{S}) where {T, S} = r + (-1 * p)

    Base.:-(p::Poly{S}, r::RationalPolynomial{T}) where {T, S} = p + (-1 * r)
    
    Base.:/(c::S, r::RationalPolynomial{T}) where {T, S<:Number} = c * inverse(r)

    Base.:/(r::RationalPolynomial{T}, c::S) where {T, S<:Number} = r * (1.0 / c)

    Base.:/(r1::RationalPolynomial{T}, r2::RationalPolynomial{T}) where {T} = r1 * inverse(r2)

    Base.:/(r::RationalPolynomial{T}, p::Poly{S}) where {T, S} = reduce(RationalPolynomial(r.numerator, r.denominator * p))

    Base.:/(p::Poly{S}, r::RationalPolynomial{T}) where {T, S} = p * inverse(r)

    Base.one(r::RationalPolynomial{T}) where {T} = RationalPolynomial(one(r.numerator), one(r.denominator))

    Base.zero(r::RationalPolynomial{T}) where {T} = RationalPolynomial(zero(r.numerator), one(r.denominator))
end
