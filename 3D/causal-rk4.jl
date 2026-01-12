"""
causal-rk4

This Julia script generates electric and magnetic fields using a causal
perspective on Maxwell's equations, and a Runge-Kutta 4 solver to estimate
solutions.

The result is four .jld2 files stored in ./cache that can be plotted
using plot.jl. There is one .jl2d file each for the charge positions,
electric fields, magnetic fields, and grid points.

contraflux
9/26/2025
"""

# ------------------------------ Dependencies ------------------------------ #

using JLD2
using LinearAlgebra

include("../calculus.jl")
include("../physics.jl")

# ------------------------------ Parameters ------------------------------ #

const dt = 0.1 # Time step 
const steps = 100 # Number of steps
const xs = collect(LinRange(-5, 5, 15))
const ys = collect(LinRange(-5, 5, 15))
const zs = collect(LinRange(-5, 5, 15))
const grid = [[x, y, z] for x in xs for y in ys for z in zs]

# ------------------------------ Point Charges ------------------------------ #

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

# ------------------------------ Simulation ------------------------------ #

function saveToBinary(grid, chargePositions, electricFields, magneticFields)
    jldsave("3D/cache/grid.jld2"; grid=grid)
    jldsave("3D/cache/charges.jld2"; charges=chargePositions)
    jldsave("3D/cache/electric.jld2"; E=electricFields)
    jldsave("3D/cache/magnetic.jld2"; B=magneticFields)
end

simulate()