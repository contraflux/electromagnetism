"""
causal-rk4

This Julia script generates electric and magnetic fields using a causal
perspective on Maxwell's equations, and a Runge-Kutta 4 solver to estimate
solutions.

The result is four .jld2 files stored in ./cache that can be plotted
using plot.jl. There is one .jl2d file each for the charge positions,
electric fields, magnetic fields, and grid points.

contraflux
1/12/2026
"""

# ------------------------------ Dependencies ------------------------------ #

using JLD2
using LinearAlgebra
using REPL
using Interpolations

include("../calculus.jl")
include("../physics.jl")

# ------------------------------ Parameters ------------------------------ #

const dt = 0.1 # Time step 
const steps = 50 # Number of steps
const xs = collect(LinRange(-5, 5, 10))
const ys = collect(LinRange(-5, 5, 10))
const zs = collect(LinRange(-5, 5, 10))
const grid = [[x, y, z] for x in xs for y in ys for z in zs]

# ------------------------------ Point Charges ------------------------------ #

function getPointCharges()
    return [PointCharge(t -> [sin(t), 0, 0], 1)]
end

# ------------------------------ Simulation ------------------------------ #

function simulate()
    pointCharges = getPointCharges()
    runtimes = []
    chargePositions = []
    electricFields = []
    magneticFields = []
    E = [zeros(3) for _ in xs for _ in ys for _ in zs]
    B = [zeros(3) for _ in xs for _ in ys for _ in zs]

    println("Running simulation with $steps steps and $dt timestep...")

    println()
    terminal = REPL.Terminals.TTYTerminal(string(), stdin, stdout, stderr)
    for t in 0:dt:(dt * steps)
        push!(runtimes, time())

        E, B = evolve(E, B, pointCharges, t)

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

function evolve(E, B, pointCharges, t)
    function electricField(x)
        e = zeros(3)
        for pc in pointCharges
            r̂(t) = normalize( x - pc.x(t) )
            r(t) = hypot( (x - pc.x(t))... )

            t′ = tprime(t -> x, pc.x, t, c)
            r̂_t′ = r̂(t′)
            r_t′ = r(t′)

            e += pc.q * r̂_t′ / (4 * π * ϵ₀ * r_t′^2)
        end
        return e
    end

    E_new = []
    B_new = deepcopy(B)

    # Static electric field
    for x in xs for y in ys for z in zs
        e = electricField([x, y, z])

        push!(E_new, e)
    end end end

    # Faraday's Law
    E_function = vectorFieldInterpolation(E)
    for i in eachindex(xs) for j in eachindex(ys) for k in eachindex(zs)
        pos = [xs[i], ys[j], zs[k]]
        e_curl = curl(E_function, pos)
        B_new[getIndex(i, j, k)] -= e_curl * dt
    end end end

    # Ampere's Law with Maxwell's Correction
    B_function = vectorFieldInterpolation(B)
    for i in eachindex(xs) for j in eachindex(ys) for k in eachindex(zs)
        pos = [xs[i], ys[j], zs[k]]
        b_curl = curl(B_function, pos)
        radius = 0.5
        J = zeros(3)
        # Approximate charge as a small sphere with constant density
        for pc in pointCharges
            ρ = (3/(4*π)) * pc.q * (1/radius^3)
            if (hypot( (pc.x(t) - pos)... ) <= radius)
                J += ρ * d_dt(pc.x, t)
            end
        end
        E_new[getIndex(i, j, k)] += (1/ϵ₀) * (b_curl - J)
    end end end
    return E_new, B_new
end

function getIndex(i, j, k)
    return i + (j-1)*length(xs) + (k-1)*length(xs)*length(ys)
end

function vectorFieldInterpolation(F)
    functions = []
    x_range = xs[begin]:xs[2]-xs[1]:xs[end]
    y_range = ys[begin]:ys[2]-ys[1]:ys[end]
    z_range = zs[begin]:zs[2]-zs[1]:zs[end]

    for i in 1:3
        itp = interpolate(reshape([f[i] for f in F], length(xs), length(ys), length(zs)), BSpline(Linear()))
        sitp = scale(itp, x_range, y_range, z_range)
        push!(functions, extrapolate(sitp, Line()))
    end

    return v -> [functions[i](v...) for i in eachindex(functions)]
end

function saveToBinary(grid, chargePositions, electricFields, magneticFields)
    jldsave("3D/cache/grid.jld2"; grid=grid)
    jldsave("3D/cache/charges.jld2"; charges=chargePositions)
    jldsave("3D/cache/electric.jld2"; E=electricFields)
    jldsave("3D/cache/magnetic.jld2"; B=magneticFields)
end

simulate()