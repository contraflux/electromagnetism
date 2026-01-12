"""
physics

This Julia script contains structs and functions used in both the causal-rk4
and jefimenk-feynman simulations.

contraflux
9/26/2025
"""

const ϵ₀ = 1
const c = 1

mutable struct PointCharge
    x::Function
    q::Float64
end

function tprime(p, x, t, c)
    δt = 0.05
    Δt = 0
    lastInLightcone = false
    scaling = 0.5

    for _ in 1:100
        distance = hypot((x(t - Δt) - p(t))...)
        lightcone = c * Δt
        nowInLightcone = distance <= lightcone
        
        if (lastInLightcone == false && nowInLightcone == true)
            δt = abs(δt) * scaling
        elseif (lastInLightcone == true && nowInLightcone == false)
            δt = -abs(δt) * scaling
        end

        Δt += δt
    end

    return t - Δt
end