%truncation_error(0)

backward_euler(@rate_func01, 0.5, solution01(0.5), 0.01)

function dXdt = rate_func01(t,X)
    dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
    X = cos(t);
end