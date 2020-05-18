using Polynomials
using LinearAlgebra

export findMax, equidistant

function equidistant(a::Number, b::Number, n::Unsigned)
    step = (b - a) / (n - 1)
    return collect(a:step:b)
end

function findMax(f::Function, a::Number, b::Number, nSteps::Unsigned = unsigned(10000), depth::Unsigned = unsigned(0))
    if a == b || depth > 1000
        return (a + b) / 2.0
    end
    step = (b - a) / (nSteps - 1)
    ys = map(x -> f(x), a:step:b)
    maxy, maxi = findmax(ys)
    maxx = a + step * (maxi - 1)
    maxx = min(b - step, max(a + step, maxx))
    findMax(f, maxx - step, maxx + step, unsigned(5), depth + 1)
end

function getAlternatingOnes(ones::AbstractArray{<:Number})
    polyfit(ones, [(-1) ^ i for i in 1:length(ones)])
end

function getPseudoAlternance(f::Function, points::AbstractArray{<:Number})
    n = length(points)
    ys = map(x -> f(x), points[1:end - 1])
    base = polyfit(points[1:end - 1], ys)
    fitting = getAlternatingOnes(points[1:end - 1])
    xn = points[end]
    e1, e2 = base(xn) - f(xn), fitting(xn)
    return base + fitting * e1 / (((-1) ^ n) - e2)
end

function getBetterFitting(f::Function, points::AbstractArray{<:Number}, a::Number, b::Number, previous::Poly, maxx::Number)
    ny, nfy = previous(maxx), f(maxx)
    if maxx <= points[1]
        py, pfy = previous(points[1]), f(points[1])
        if sign(py - pfy) != sign(ny - nfy)
            points = circshift(points, 1)
        end
        points[1] = maxx
        return points
    elseif maxx >= points[end]
        py, pfy = previous(points[end]), f(points[end])
        if sign(py - pfy) != sign(ny - nfy)
            points = circshift(points, -1)
        end
        points[end] = maxx
        return points
    else
        left = findfirst(x -> maxx <= x, points) - 1
        py, pfy = previous(points[left]), f(points[left])
        if sign(py - pfy) == sign(ny - nfy)
            points[left] = maxx
        else
            points[left + 1] = maxx
        end
        return points
    end
end

function getAlternance(f::Function, from::Real, to::Real, degree::Unsigned, precision::Real = 10e-20, limit::Unsigned = unsigned(100))
    if precision <= 0
        throw(DomainError(precision, "precision must be greater than 0"))
    elseif from >= to
        throw(Exception("to has to be greater than from"))
    end
    points = equidistant(from, to, degree + 2)
    alternance = getPseudoAlternance(f, points)
    maxx = findMax(x -> abs(f(x) - alternance(x)), from, to)
    error = abs(f(maxx) - alternance(maxx)) - abs(f(points[1]) - alternance(points[1]))
    iter = 0
    while error >= precision || iter >= limit
        iter += 1
        points = getBetterFitting(f, points, from, to, alternance, maxx)
        alternance = getPseudoAlternance(f, points)
        maxx = findMax(x -> abs(f(x) - alternance(x)), from, to)
        error = abs(f(maxx) - alternance(maxx)) - abs(f(points[1]) - alternance(points[1]))
    end
    return points, alternance
end


f = x -> sin(x)
a = 0.0
b = 1.0
s = x -> x - ((x^3)/6) + ((x^5)/(120))
points, fa = getAlternance(f, a, b, unsigned(3))
using Plots
    gr()
    points = a:((b - a) / 1000):b
    plot(x -> abs(fa(x) - f(x)), points)
    plot!(x -> abs(s(x) - f(x)), points)

function primitiveTest()
    return ["using Plots
    gr()
    f = x -> sin(3x)
    points = [i for i in -2.0:0.5:2.0]

    fa = getAlternance(f, points)
    plot(fa, -2.0:0.001:2.0)
    plot!(f, -2.0:0.001:2.0)
    scatter!(fa, points)
    points = getBetterFitting(f, points, -2.0, 2.0)"]
end


