const δ = 1e-5

function d_dt(f, t)
    ( f(t + δ) - f(t) ) / δ
end

function ∂_∂x(F, v)
    ∂x = [δ, 0, 0]
    ( F(v + ∂x) - F(v) ) / δ
end

function ∂_∂y(F, v)
    ∂y = [0, δ, 0]
    ( F(v + ∂y) - F(v) ) / δ
end

function ∂_∂z(F, v)
    ∂z = [0, 0, δ]
    ( F(v + ∂z) - F(v) ) / δ
end

function div(F, v)
    ∂_∂x(F, v)[1] + ∂_∂y(F, v)[2] + ∂_∂z(F, v)[3]
end

function curl(F, v)
    i = ∂_∂y(F, v)[3] - ∂_∂z(F, v)[2]
    j = ∂_∂z(F, v)[1] - ∂_∂x(F, v)[3]
    k = ∂_∂x(F, v)[2] - ∂_∂y(F, v)[1]
    [i, j, k]
end