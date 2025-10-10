%truncation_error(0)

[t_list, x_list, h_avg, num_evals] = backward_euler(@rate_func01, [0,(pi/4)], 1, 0.01)
plot(t_list,x_list); hold on
plot(t_list, solution01(t_list))

function dXdt = rate_func01(t,X)
    dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
    X = cos(t);
end