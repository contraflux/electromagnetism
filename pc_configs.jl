function getPointCharges()
    function toroid(t, ϕ)
        w = 0.5
        r = 2.5
        k = 16

        s = 0.075*c*t + ϕ

        return [cos(s) * (r + w * sin(k * s)), sin(s) * (r + w * sin(k * s)), w * cos(k * s)]
    end

    return [PointCharge(t -> toroid(t, ϕ), 1) for ϕ in 0:π/16:31π/16]
end