f(E, M, e) = E - e * sin(E) - M


df(E, e) = 1 - e * cos(E)

newton_rapshon_iteration(E, M, e) = E - f(E, M, e) / df(E, e)

function newton_raphson(M, e; threshold = 1e-5)
    E = M
    i = 0
    loss = f(E, M, e)

    while abs(loss) > threshold && i < 100
        E = newton_rapshon_iteration(E, M, e)
        i += 1
        loss = f(E, M, e)
    end

    return E, loss
end
