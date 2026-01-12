using GLMakie

function rk4_step(ẋ, δt, xₙ, tₙ)
    k1 = ẋ(tₙ, xₙ)
    k2 = ẋ(tₙ + δt/2, xₙ + (δt * k1/2))
    k3 = ẋ(tₙ + δt/2, xₙ + (δt * k2/2))
    k4 = ẋ(tₙ + δt, xₙ + (δt * k3))

    xₙ₊₁ = xₙ + (δt / 6) * (k1 + 2k2 + 2k3 + k4)
    return xₙ₊₁
end

function runge_kutta(ẋ, δt, x₀, t₀, iterations)
    xₙ = x₀
    tₙ = t₀

    xs = [];

    for n in 1:iterations
        append!(xs, xₙ)
        xₙ = rk4_step(ẋ, δt, xₙ, tₙ)
        tₙ = tₙ + δt
    end

    return xs
end

ẋ(t, x) = 2t^2 - 4t + t

n = 100

ts = collect(LinRange(0, 3, n))
δt = ts[2] - ts[1]
xs = runge_kutta(ẋ, δt, 0, 0, n)

fig = Figure()
ax = Axis(fig[1, 1])

lines!(ax, ts, xs, color=:black)
display(fig)