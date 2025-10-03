%Runs numerical integration using forward Euler approximation
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
[t_list,X_list,h_avg, num_evals] = forward_euler_test(@rate_func01, [0, pi], 0, pi/4);
X_list = X_list
X = solution01(t_list)

function [t_list,X_list,h_avg, num_evals] = forward_euler_test(rate_func_in,tspan,X0,h_ref)
    num_of_steps = ceil((tspan(2) - tspan(1)) / h_ref);
    h_avg = (tspan(2) - tspan(1)) / num_of_steps;
    t_list = linspace(tspan(1), tspan(2), num_of_steps)';
    X_list = [X0];
    num_evals = 0;

    for i = 2:length(t_list)
        [XB_temp,num_eval_temp] = forward_euler_step(rate_func_in, t_list(i), X_list(i-1), h_avg);
        X_list = [X_list; XB_temp];
        num_evals = num_evals + num_eval_temp;
    end

end


function [XB,num_evals] = forward_euler_step(rate_func_in,t,XA,h)
    % computing for value of x when t increases by h
    XB = XA + h*rate_func_in(t, XA);

    % account for how many times we called rate_func_in
    num_evals = 1;
end

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