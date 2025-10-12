

t0 = pi/6;
h = logspace(-5, 1, 100);
analytical_difference = zeros(1, length(h));
exp_euler_local_error = zeros(1, length(h));
exp_midpoint_local_error = zeros(1, length(h));
imp_euler_local_error = zeros(1, length(h));
imp_midpoint_local_error = zeros(1, length(h));

h_avg = zeros(1, length(h));

for i = 1:length(h)
    analytical_difference(i) = norm(solution01(t0+h(1,i))-solution01(t0));
    [exp_euler_local_error(i), ~, ~] = forward_euler_local_error(t0, h(1,i));
    [exp_midpoint_local_error(i), h_avg_temp, ~] = explicit_midpoint_local_error(t0, h(1,i));
    [imp_euler_local_error(i), ~, ~] = backward_euler_local_error(t0, h(1,i));
    [imp_midpoint_local_error(i), ~, ~] = implicit_midpoint_local_error(t0, h(1,i));

    h_avg(1, i) = h_avg_temp;
end

[p1,k1] = loglog_fit(h,analytical_difference);
[p2,k2] = loglog_fit(h_avg,exp_euler_local_error);
[p3,k3] = loglog_fit(h_avg,exp_midpoint_local_error);
[p4,k4] = loglog_fit(h_avg,imp_euler_local_error);
[p5,k5] = loglog_fit(h_avg,imp_midpoint_local_error);
disp(p1)
disp(p2)
disp(p3)
disp(p4)
disp(p5)


function [local_error, h_avg, num_evals] = forward_euler_local_error(t, h)
    % calculate local error for the first test function
    [~,X_list,h_avg, num_evals] = forward_euler(@rate_func01,[t, t+h],solution01(t),h+1);
    G_t = X_list(end);
    X_t = solution01(t+h);

    local_error = norm(G_t - X_t);
end

function [local_error, h_avg, num_evals] = explicit_midpoint_local_error(t, h)
    % calculate local error for the first test function
    
    [~,X_list, h_avg, num_evals] = explicit_midpoint(@rate_func01,[t, t+h],solution01(t),h+1);
    
    G_t = X_list(end);
    X_t = solution01(t+h);

    local_error = norm(G_t - X_t);
end

function [local_error, h_avg, num_evals] = backward_euler_local_error(t, h)
    % calculate local error for the first test function
    [~,X_list,h_avg, num_evals] = backward_euler(@rate_func01,[t, t+h],solution01(t),h+1);
    G_t = X_list(end);
    X_t = solution01(t+h);

    local_error = norm(G_t - X_t);
end

function [local_error, h_avg, num_evals] = implicit_midpoint_local_error(t, h)
    % calculate local error for the first test function
    
    [~,X_list, h_avg, num_evals] = implicit_midpoint(@rate_func01,[t, t+h],solution01(t),h+1);
    
    G_t = X_list(end);
    X_t = solution01(t+h);

    local_error = norm(G_t - X_t);
end

function [global_error, h_avg, num_evals] = forward_euler_global_error(tspan, h_ref)
    % calculate global error for the first test function
    [t_list,X_f,h_avg, num_evals] = forward_euler(@rate_func01,tspan,solution01(tspan(1)),h_ref);
    X_tf = solution01(t_list)';

    global_error = norm(X_f - X_tf);
end

function [global_error, h_avg, num_evals] = explicit_midpoint_global_error(tspan, h_ref)
    % calculate global error for the first test function
    [t_list,X_f, h_avg, num_evals] = explicit_midpoint(@rate_func01,tspan,solution01(tspan(1)),h_ref);
    X_tf = solution01(t_list)';

    global_error = norm(X_f - X_tf);
end

function [global_error, h_avg, num_evals] = backward_euler_global_error(tspan, h_ref)
    % calculate global error for the first test function
    [t_list,X_f,h_avg, num_evals] = backward_euler(@rate_func01,tspan,solution01(tspan(1)),h_ref);
    X_tf = solution01(t_list)';

    global_error = norm(X_f - X_tf);
end

function [global_error, h_avg, num_evals] = implicit_midpoint_global_error(tspan, h_ref)
    % calculate global error for the first test function
    [t_list,X_f, h_avg, num_evals] = implicit_midpoint(@rate_func01,tspan,solution01(tspan(1)),h_ref);
    X_tf = solution01(t_list)';

    global_error = norm(X_f - X_tf);
end

% test funcs_______________________________________________________________
% test func 1
% function dXdt = rate_func01(t,X)
%     dXdt = -5*X + 5*cos(t) - sin(t);
% end
% 
% function X = solution01(t)
%     X = cos(t);
% end

% test func 2
function dXdt = rate_func01(t,X)
    dXdt = [0,-1;1,0]*X;
end

function X = solution01(t)
    X = [cos(t);sin(t)];
end