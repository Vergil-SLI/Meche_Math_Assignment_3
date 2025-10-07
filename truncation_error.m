function truncation_error()
    clear;
    % % plots comparing closed form approximation (of various time step) with numerical solutions
    % % forward euler method
    h = logspace(0.01, 0.1, 5)
    t_start = 0;
    t_end = pi;

    [t_list, euler_X_list, h_avg, ~] = forward_euler(@rate_func01,[t_start, t_end],solution01(t_start), h(1));
    
    axis equal
    hold off
    numerical_X_list = solution01(t_list);
    plot(t_list, numerical_X_list, 'black', LineWidth=2)

    hold on
    for i = 1:length(h)
        [t_list, euler_X_list, h_avg, ~] = forward_euler(@rate_func01,[0, pi],solution01(t_start), h(i));
        
        plot(t_list, euler_X_list, "o-", LineWidth=2)
    end
    
    legend()
    

    % plot comparing the 4 different approximation methods when given the same time step
    


    % h = logspace(-5, 1, 100);
    % euler_local_error = zeros(1, length(h));
    % midpoint_local_error = zeros(1, length(h));
    % 
    % for i = 1:length(h)
    %     euler_local_error(i) = forward_euler_local_error(pi/2, h(1,i));
    %     midpoint_local_error(i) = explicit_midpoint_local_error(pi/2, h(1,i));
    % end
    % 
    % % plot local truncation error
    % axis equal
    % hold off
    % loglog(h,euler_local_error, 'b-')
    % hold on
    % loglog(h,midpoint_local_error, 'g-')
    % legend("euler forward local error", "explicit midpoint local error")

    % 
    % loglog_fit(h,euler_local_error)
    % loglog_fit(h,midpoint_local_error)
    % 
    % % global_error = forward_euler_global_error([0, 1], 0.01);
    % % global_error = explicit_midpoint_global_error([0, 1], 0.01);
end


function local_error = forward_euler_local_error(t, h)
    % calculate local error for the first test function
    [~,X_list,~, ~] = forward_euler(@rate_func01,[t, t+h],solution01(t),h+1);
    G_t = X_list(end);
    X_t = solution01(t+h);

    local_error = norm(G_t - X_t);
end

function local_error = explicit_midpoint_local_error(t, h)
    % calculate local error for the first test function
    [~,X_list,~, ~] = explicit_midpoint(@rate_func01,[t, t+h],solution01(t),h+1);
    G_t = X_list(end);
    X_t = solution01(t+h);

    local_error = norm(G_t - X_t);
end

function global_error = forward_euler_global_error(tspan, h_ref)
    % calculate global error for the first test function
    [t_list,X_f,~, num_evals] = forward_euler(@rate_func01,tspan,solution01(tspan(1)),h_ref);
    X_tf = solution01(t_list)';

    global_error = norm(X_f - X_tf);
end

function global_error = explicit_midpoint_global_error(tspan, h_ref)
    % calculate global error for the first test function
    [t_list,X_f,~, num_evals] = explicit_midpoint(@rate_func01,tspan,solution01(tspan(1)),h_ref);
    X_tf = solution01(t_list)';

    global_error = norm(X_f - X_tf);
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