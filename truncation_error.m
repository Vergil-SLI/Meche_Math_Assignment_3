function truncation_error()
    forward_euler_local_error(2, 0.001)
    explicit_midpoint_local_error(2, 0.001)
end

function local_error = forward_euler_local_error(t, h)
    G_t = forward_euler_step(@rate_func01,t,solution01(t),h);
    X_t = solution01(t+h);

    local_error = norm(G_t - X_t);
end

function local_error = explicit_midpoint_local_error(t, h)
    G_t = explicit_midpoint_step(@rate_func01,t,solution01(t),h);
    X_t = solution01(t+h);

    local_error = norm(G_t - X_t);
end

% explicit midpoint________________________________________________________
function [t_list,X_list,h_avg, num_evals] = explicit_midpoint(rate_func_in,tspan,X0,h_ref)
    num_of_steps = ceil((tspan(2) - tspan(1)) / h_ref);
    h_avg = (tspan(2) - tspan(1)) / num_of_steps;
    t_list = linspace(tspan(1), tspan(2), num_of_steps+1)';
    X_list = zeros(length(X0),length(t_list));
    X_list(:,1) = X0;
    num_evals = 0;

    for i = 2:length(t_list)
        [XB_temp,num_eval_temp] = explicit_midpoint_step(rate_func_in, t_list(i-1), X_list(i-1), h_avg);
        X_list(:,i) = XB_temp;
        num_evals = num_evals + num_eval_temp;
    end
end

%This function computes the value of X at the next time step
%using the explicit midpoint approximation
function [XB,num_evals] = explicit_midpoint_step(rate_func_in,t,XA,h)
    XB_half = XA + (h/2) * rate_func_in(t, XA);
    XB = XB_half + (h/2) * rate_func_in(t+(h/2), XB_half);
    
    num_evals = 2;
end


% forward_euler____________________________________________________________
function [t_list,X_list,h_avg, num_evals] = forward_euler(rate_func_in,tspan,X0,h_ref)
    num_of_steps = ceil((tspan(2) - tspan(1)) / h_ref);
    h_avg = (tspan(2) - tspan(1)) / num_of_steps;
    t_list = linspace(tspan(1), tspan(2), num_of_steps+1)';
    X_list = zeros(length(X0),length(t_list));
    X_list(:,1) = X0;
    num_evals = 0;

    for i = 2:length(t_list)
        [XB_temp,num_eval_temp] = forward_euler_step(rate_func_in, t_list(i-1), X_list(i-1), h_avg);
        X_list(:,i) = XB_temp;
        num_evals = num_evals + num_eval_temp;
    end

end

function [XB,num_evals] = forward_euler_step(rate_func_in,t,XA,h)
    % computing for value of x when t increases by h
    XB = XA + h*rate_func_in(t, XA);

    % account for how many times we called rate_func_in
    num_evals = 1;
end


% test funcs_______________________________________________________________
% test func 1
function dXdt = rate_func01(t,X)
    dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
    X = cos(t);
end

% test func 2
function dXdt = rate_func02(t,X)
    dXdt = [0,-1;1,0]*X;
end

function X = solution02(t)
    X = [cos(t);sin(t)];
end