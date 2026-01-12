"""
jefimenko-feynman

This Julia script generates electric and magnetic fields using the
Jefimenko-Feynman (also known as Heaviside-Feynman) equations.

The result is four .jld2 files stored in ./cache that can be plotted
using plot.jl. There is one .jl2d file each for the charge positions,
electric fields, magnetic fields, and grid points.

contraflux
9/26/2025
"""

# ------------------------------ Dependencies ------------------------------ #

using JLD2
using LinearAlgebra
using REPL

include("../calculus.jl")
include("../physics.jl")

# ------------------------------ Parameters ------------------------------ #

const dt = 0.1 # Time step 
const steps = 300 # Number of steps
const xs = collect(LinRange(-5, 5, 15))
const ys = collect(LinRange(-5, 5, 15))
const zs = collect(LinRange(-5, 5, 15))
const grid = [[x, y, z] for x in xs for y in ys for z in zs]

# ------------------------------ Point Charges ------------------------------ #

function getPointCharges()
    return [PointCharge(t -> [1.5*q*cos(c*0.6*t), 0, 0], q) for q in -1:2:1]
end

# ------------------------------ Simulation ------------------------------ #

function simulate()
    pointCharges = getPointCharges()
    runtimes = []
    chargePositions = []
    electricFields = []
    magneticFields = []

    println("Running simulation with $steps steps and $dt timestep...")

    println()
    terminal = REPL.Terminals.TTYTerminal(string(), stdin, stdout, stderr)
    for t in 0:dt:(dt * steps)
        push!(runtimes, time())

        E, B = jefimenkoFeynman(pointCharges, t)

        push!(electricFields, E)
        push!(magneticFields, B)
        push!(chargePositions, [pc.x(t) for pc in pointCharges])

        text = "$(round(100 * t / (dt * steps), digits=3))%"
        REPL.Terminals.cmove_left(terminal, length(text))
        REPL.Terminals.cmove_line_up(terminal)
        println(text)
    end

    total_runtime = round(runtimes[end] - runtimes[begin], digits=3)
    runtime_per_step = round(total_runtime / steps, digits=3)
    println("Run completed in $total_runtime seconds, with $runtime_per_step seconds per step.")

    println("Saving run to the cache...")
    saveToBinary(grid, chargePositions, electricFields, magneticFields)
    println("Saved to the cache!")
end

function jefimenkoFeynman(pointCharges, t)
    function electromagneticField(x)
        e = zeros(3)
        b = zeros(3)
    
        for pc in pointCharges
            r̂(t) = normalize( x - pc.x(t) )
            r(t) = hypot((  x - pc.x(t) )... )

            t′ = tprime(t -> x, pc.x, t, c)
            r̂_t′ = r̂(t′)
            r_t′ = r(t′)

            e_1 = r̂_t′ / r_t′^2
            e_2 = r_t′ * d_dt(t -> (r̂(t) / r(t)^2), t′) / c
            e_3 = d_dt(t -> d_dt(r̂, t), t′) / c^2
            e_q = -pc.q * (e_1 + e_2 + e_3) / (4 * π * ϵ₀)

            e += e_q 
            b += cross(-r̂_t′, e_q / c)
        end

        return e, b
    end

    E = []
    B = []

    for x in xs for y in ys for z in zs
        e, b = electromagneticField([x, y, z])

        push!(E, e)
        push!(B, b)
    end end end

    return E, B
end

function saveToBinary(grid, chargePositions, electricFields, magneticFields)
    jldsave("3D/cache/grid.jld2"; grid=grid)
    jldsave("3D/cache/charges.jld2"; charges=chargePositions)
    jldsave("3D/cache/electric.jld2"; E=electricFields)
    jldsave("3D/cache/magnetic.jld2"; B=magneticFields)
end

simulate()