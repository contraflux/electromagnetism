"""
plot

This Julia script generates a plot and .mp4 video of the electric and magnetic
fields generated using either the causal-rk4 or jefimenk-feynman scripts.

contraflux
9/26/2025
"""

# ------------------------------ Dependencies ------------------------------ #

using GLMakie
using JLD2
using REPL

# ------------------------------ Parameters ------------------------------ #

const framerate = 30
const fig = Figure(size=(800, 600))
const ax = Axis3(fig[1, 1], limits=(-5, 5, -5, 5, -5, 5))

# ------------------------------ Plotting ------------------------------ #

function plot(filename)
    set_theme!(merge(theme_black(), theme_latexfonts()))

    println("Loading run from the cache...")
    grid, charges, Es, Bs = loadFromBinary()
    println("Loaded from the cache!")

    points = [Point3f(g) for g in grid]

    println("Plotting simulation...")

    println()
    terminal = REPL.Terminals.TTYTerminal(string(), stdin, stdout, stderr)
    record(fig, filename * ".mp4", 1:length(Es); framerate = framerate) do i
        empty!(ax)

        E = Es[i]
        B = Bs[i]
        charge_positions = charges[i]

        electric_lengths = [hypot(e...) for e in E]
        electric_max_length = maximum(electric_lengths)
        electric_colors = [RGBAf(1, 0, 0, l) for l in electric_lengths]

        magnetic_lengths = [hypot(b...) for b in B]
        magnetic_max_length = maximum(magnetic_lengths)
        magnetic_colors = [RGBAf(0, 1e1 * l, 2e1 * l, 1e1 * l) for l in magnetic_lengths]

        arrows!(ax, points, [Vec3f(e / electric_max_length) for e in E], normalize=true, color=electric_colors, lengthscale=0.1, arrowsize=0.015, transparency=true)
        arrows!(ax, points, [Vec3f(b / magnetic_max_length) for b in B], normalize=true, color=magnetic_colors, lengthscale=0.1, arrowsize=0.015, transparency=true)

        for position in charge_positions
            scatter!(ax, position..., color=:white, markersize=8)
        end

        display(fig)

        if (i % 10 == 0)
            text = "Frame $i of $(length(Es))"
            REPL.Terminals.cmove_left(terminal, length(text))
            REPL.Terminals.cmove_line_up(terminal)
            println(text)
        end
        
        ax.elevation[] = π/8
        ax.azimuth[] = (π/8) + (π/4) * i / (101 * 2)
    end

    println("Plot complete! Saved to \" " * filename * ".mp4\" ")
end

function loadFromBinary()
    grid = load("3D/cache/grid.jld2", "grid")
    charges = load("3D/cache/charges.jld2", "charges")
    Es = load("3D/cache/electric.jld2", "E")
    Bs = load("3D/cache/magnetic.jld2", "B")

    return grid, charges, Es, Bs
end

plot("sim_run")