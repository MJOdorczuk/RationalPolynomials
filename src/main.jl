
include("RationalFunctionAnalysis.jl")
include("PolynomialUtils.jl")

p = 2 * poly([1., 2.])
q = 4 * poly([1., 2., 3.])
r = RationalPolynomials.RationalPolynomial(p, q)
PolyUtils.polyLongDivision(p, q)