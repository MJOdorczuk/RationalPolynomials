module PolyUtils
    export polyLongDivisionWithReminder, polyLongDivision
    using Polynomials
    polyX = Poly([0, 1])

    function polyLongDivisionWithReminder(p1::Poly{T}, p2::Poly{S}) where {T, S}
        deg1, deg2 = degree(p1), degree(p2)
        if deg1 < deg2
            return (Poly([0]), p1)
        elseif deg1 == deg2
            a = p1.a[end] / p2.a[end]
            q = Poly([a])
            return (q, p1 - a * p2)
        else
            q, r = polyLongDivisionWithReminder(p1, p2 * polyX)
            a = r.a[end] / p2.a[end]
            q = q * polyX + a
            return (q, r - a * p2)
        end
    end

    function polyLongDivision(p1::Poly{T}, p2::Poly{S}) where {T, S} 
        return polyLongDivisionWithReminder(p1, p2)[1]
    end

end
