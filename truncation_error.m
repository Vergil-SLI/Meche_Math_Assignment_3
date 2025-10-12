function truncation_error(mode)
    if mode == 1
        % Comparing closed form approximation (of various time step) with numerical solutions
        % Forward Euler Method
        h = linspace(0.375, 0.475, 3); % WHAT ARE GOOD H VALUES
        t_start = 0;
        t_end = 7*pi/4;
        legend_titles = cell(1, length(h));
        legend_titles{1, 1} = "numerical solution";

        % plot the numerical solution
        hold off
        t_list = linspace(t_start, t_end, 20);
        analytical_X_list = solution01(t_list);
        plot(t_list, analytical_X_list, 'black')
        hold on

        % plot the approximations (of various time step)
        for i = 1:length(h)
            [t_list, foward_euler_X_list, h_avg, ~] = forward_euler(@rate_func01,[t_start, t_end],solution01(t_start), h(i));
            plot(t_list, foward_euler_X_list, "-")
            legend_titles{1, i+1} = "closed form approx (h = " + num2str(round(h_avg, 3)) + ")";
        end
        lgd = legend(legend_titles);
        lgd.Location = "southeast";
        xlabel("t")
        ylabel("X Value")
        title("Comparing Forward Euler Approximation with Analytical Solution")
    end

    % _____________________________________________________________________
    % comparing explicit midpoint with closed form solution
    if mode == 2
       % Comparing closed form approximation (of various time step) with numerical solutions
        % Explicit Midpoint
        h = linspace(0.375, 0.475, 3); % WHAT ARE GOOD H VALUES
        t_start = 0;
        t_end = 7*pi/4;
        legend_titles = cell(1, length(h));
        legend_titles{1, 1} = "numerical solution";

        % plot the numerical solution
        hold off
        t_list = linspace(t_start, t_end, 20);
        analytical_X_list = solution01(t_list);
        plot(t_list, analytical_X_list, 'black')
        hold on

        % plot the approximations (of various time step)
        for i = 1:length(h)
            [t_list, foward_euler_X_list, h_avg, ~] = explicit_midpoint(@rate_func01,[t_start, t_end],solution01(t_start), h(i));
            plot(t_list, foward_euler_X_list, "-")
            legend_titles{1, i+1} = "closed form approx (h = " + num2str(round(h_avg, 3)) + ")";
        end
        lgd = legend(legend_titles);
        lgd.Location = "southeast";
        xlabel("t")
        ylabel("X Value")
        title("Comparing Explicit Midpoint Approximation with Analytical Solution")
    end

    % _____________________________________________________________________
    % comparing implicit euler with closed form solution
    if mode == 3
       % Comparing closed form approximation (of various time step) with numerical solutions
        % Backward Euler Method
        h = linspace(0.375, 0.475, 3); % WHAT ARE GOOD H VALUES
        t_start = 0;
        t_end = 7*pi/4;
        legend_titles = cell(1, length(h));
        legend_titles{1, 1} = "numerical solution";

        % plot the numerical solution
        hold off
        t_list = linspace(t_start, t_end, 20);
        analytical_X_list = solution01(t_list);
        plot(t_list, analytical_X_list, 'black')
        hold on

        % plot the approximations (of various time step)
        for i = 1:length(h)
            [t_list, foward_euler_X_list, h_avg, ~] = backward_euler(@rate_func01,[t_start, t_end],solution01(t_start), h(i));
            plot(t_list, foward_euler_X_list, "-")
            legend_titles{1, i+1} = "closed form approx (h = " + num2str(round(h_avg, 3)) + ")";
        end
        lgd = legend(legend_titles);
        lgd.Location = "southeast";
        xlabel("t")
        ylabel("X Value")
        title("Comparing Implicit (Backwards) Euler Approximation with Analytical Solution")
    end

    % _____________________________________________________________________
    %comparing implicit midpoint with closed form solution
    if mode == 4
       % Comparing closed form approximation (of various time step) with numerical solutions
        % Implicit Midpoint
        h = linspace(0.375, 0.475, 3); % WHAT ARE GOOD H VALUES
        t_start = 0;
        t_end = 7*pi/4;
        legend_titles = cell(1, length(h));
        legend_titles{1, 1} = "numerical solution";

        % plot the numerical solution
        hold off
        t_list = linspace(t_start, t_end, 20);
        analytical_X_list = solution01(t_list);
        plot(t_list, analytical_X_list, 'black')
        hold on

        % plot the approximations (of various time step)
        for i = 1:length(h)
            [t_list, foward_euler_X_list, h_avg, ~] = implicit_midpoint(@rate_func01,[t_start, t_end],solution01(t_start), h(i));
            plot(t_list, foward_euler_X_list, "-")
            legend_titles{1, i+1} = "closed form approx (h = " + num2str(round(h_avg, 3)) + ")";
        end
        lgd = legend(legend_titles);
        lgd.Location = "southeast";
        xlabel("t")
        ylabel("X Value")
        title("Comparing Implicit Midpoint Approximation with Analytical Solution")
    end

    % _____________________________________________________________________
    % Comparing the 4 different approximation methods when href = 0.38
    if mode == 5
        t_start = 0;
        t_end = 7*pi/4;
        h = 0.38;
        [t_list, foward_euler_X_list, ~, ~] = forward_euler(@rate_func01,[t_start, t_end],solution01(t_start), h);
        [~, explicit_midpoint_X_list, ~, ~] = explicit_midpoint(@rate_func01,[t_start, t_end],solution01(t_start), h);
        [~, backward_euler_X_list, ~, ~] = backward_euler(@rate_func01,[t_start, t_end],solution01(t_start), h);
        [~, implicit_midpoint_X_list, ~, ~] = implicit_midpoint(@rate_func01,[t_start, t_end],solution01(t_start), h);
        analytical_X_list = solution01(t_list);

        hold off
        plot(t_list, foward_euler_X_list, 'b-')
        hold on
        plot(t_list,explicit_midpoint_X_list, 'r-')
        plot(t_list,backward_euler_X_list, 'y-')
        plot(t_list,implicit_midpoint_X_list, 'm-')
        plot(t_list, analytical_X_list, 'g-')

        lgd = legend('forward euler', 'explicit midpoint', 'backwards euler', 'implicit midpoint', 'analytical solution');
        lgd.Location = "southeast";
        xlabel("t")
        ylabel("X Value")
        title("Comparison of all 4 approximation methods, href = 0.38")
    end


    % _____________________________________________________________________
    if mode == 6
        % Local truncation error between the 2 explicit methods
        t0 = pi/6;
        h = logspace(-5, 1, 100);
        analytical_difference = zeros(1, length(h));
        imp_euler_local_error = zeros(1, length(h));
        exp_midpoint_local_error = zeros(1, length(h));
        h_avg = zeros(1, length(h));

        for i = 1:length(h)
            analytical_difference(i) = abs(solution01(t0+h(1,i))-solution01(t0));
            [imp_euler_local_error(i), ~, ~] = forward_euler_local_error(t0, h(1,i));
            [exp_midpoint_local_error(i), h_avg_temp, ~] = explicit_midpoint_local_error(t0, h(1,i));

            h_avg(1, i) = h_avg_temp;
        end

        [p1,k1] = loglog_fit(h,analytical_difference);
        [p2,k2] = loglog_fit(h_avg,imp_euler_local_error);
        [p3,k3] = loglog_fit(h_avg,exp_midpoint_local_error);

        % plot the local truncation error
        hold off
        loglog(h,analytical_difference, 'b.', MarkerSize=10)
        hold on
        loglog(h_avg,imp_euler_local_error, 'r.', MarkerSize=10)
        loglog(h_avg,exp_midpoint_local_error, 'g.', MarkerSize=10)
        loglog(h, k1*(h.^p1), 'b-')
        loglog(h_avg, k2*(h_avg.^p2), 'r-')
        loglog(h_avg, k3*(h_avg.^p3), 'g-')
        xlabel("Average timestep length h")
        ylabel("Error")
        lgd = legend("analytical difference","forward euler local error", "explicit midpoint local error",...
            "k = " + k1 + ", p = " + p1,...
            "k = " + k2 + ", p = " + p2,...
            "k = " + k3 + ", p = " + p3);
        lgd.Location = "southeast";
        title("Local Truncation Error of Explicit Methods")
    end
    
    % local truncation of all 4 methods
    if mode == 7
        % Local truncation error between the 2 explicit methods
        t0 = pi/6;
        h = logspace(-5, 1, 100);
        %analytical_difference = zeros(1, length(h));
        exp_euler_local_error = zeros(1, length(h));
        exp_midpoint_local_error = zeros(1, length(h));
        imp_euler_local_error = zeros(1, length(h));
        imp_midpoint_local_error = zeros(1, length(h));
        h_avg = zeros(1, length(h));

        for i = 1:length(h)
            %analytical_difference(i) = abs(solution01(t0+h(1,i))-solution01(t0));
            [exp_euler_local_error(i), ~, ~] = forward_euler_local_error(t0, h(1,i));
            [exp_midpoint_local_error(i), h_avg_temp, ~] = explicit_midpoint_local_error(t0, h(1,i));
            [imp_euler_local_error(i), ~, ~] = backward_euler_local_error(t0, h(1,i));
            [imp_midpoint_local_error(i), ~, ~] = implicit_midpoint_local_error(t0, h(1,i));


            h_avg(1, i) = h_avg_temp;
        end

        %[p1,k1] = loglog_fit(h,analytical_difference);
        [p2,k2] = loglog_fit(h_avg,exp_euler_local_error);
        [p3,k3] = loglog_fit(h_avg,exp_midpoint_local_error);
        [p4,k4] = loglog_fit(h_avg,imp_euler_local_error);
        [p5,k5] = loglog_fit(h_avg,imp_midpoint_local_error);

        % plot the local truncation error
        hold off
        %loglog(h,analytical_difference, 'b.', MarkerSize=10)
        loglog(h_avg,exp_euler_local_error, 'r.', MarkerSize=10)
        hold on
        loglog(h_avg,exp_midpoint_local_error, 'g.', MarkerSize=10)
        loglog(h_avg,imp_euler_local_error, 'c.', MarkerSize=10)
        loglog(h_avg,imp_midpoint_local_error, 'm.', MarkerSize=10)
        %loglog(h, k1*(h.^p1), 'b-')
        %loglog(h_avg, k2*(h_avg.^p2), 'r-')
        %loglog(h_avg, k3*(h_avg.^p3), 'g-')
        %loglog(h_avg, k4*(h_avg.^p4), 'c-')
        %loglog(h_avg, k5*(h_avg.^p5), 'm-')
        xlabel("Average timestep length h")
        ylabel("Error")
        lgd = legend("forward euler local error", "explicit midpoint local error", "backward euler local error", "implicit midpoint local error");
        lgd.Location = "southeast";
        title("Local Truncation Error of All Methods")
    end

    % % _____________________________________________________________________
    % Global truncation error between the explicit methods 
    if mode == 8
        tspan = [0, pi/6];
        h = logspace(-5, 1, 100);
        imp_euler_global_error = zeros(1, length(h));
        exp_midpoint_global_error = zeros(1, length(h));
        h_avg = zeros(1, length(h));

        for i = 1:length(h)
            [imp_euler_global_error(i), ~, ~] = forward_euler_global_error(tspan, h(1,i));
            [exp_midpoint_global_error(i), h_avg_temp, ~] = explicit_midpoint_global_error(tspan, h(1,i));
            h_avg(1, i) = h_avg_temp;
        end

        [p1,k1] = loglog_fit(h_avg,imp_euler_global_error);
        [p2,k2] = loglog_fit(h_avg,exp_midpoint_global_error);

        % plot the local truncation error
        hold off
        loglog(h_avg,imp_euler_global_error, 'r.', MarkerSize=10)
        hold on
        loglog(h_avg,exp_midpoint_global_error, 'g.', MarkerSize=10)
        loglog(h_avg, k1*(h_avg.^p1), 'r-')
        loglog(h_avg, k2*(h_avg.^p2), 'g-')
        lgd = legend("forward euler global error", "explicit midpoint global error",...
            "k = " + k1 + ", p = " + p1,...
            "k = " + k2 + ", p = " + p2);
        lgd.Location = "southeast";
        xlabel("Average timestep length h")
        ylabel("Error")
        title("Global Truncation Error of Explicit Methods")
    end

    
    % % _____________________________________________________________________
    % Global truncation error between all methods 
    if mode == 9
        % Local truncation error between the 2 explicit methods
        tspan = [0, pi/6];
        h = logspace(-5, 1, 100);
        %analytical_difference = zeros(1, length(h));
        exp_euler_global_error = zeros(1, length(h));
        exp_midpoint_global_error = zeros(1, length(h));
        imp_euler_global_error = zeros(1, length(h));
        imp_midpoint_global_error = zeros(1, length(h));
        h_avg = zeros(1, length(h));

        for i = 1:length(h)
            %analytical_difference(i) = abs(solution01(t0+h(1,i))-solution01(t0));
            [exp_euler_global_error(i), ~, ~] = forward_euler_global_error(tspan, h(1,i));
            [exp_midpoint_global_error(i), h_avg_temp, ~] = explicit_midpoint_global_error(tspan, h(1,i));
            [imp_euler_global_error(i), ~, ~] = backward_euler_global_error(tspan, h(1,i));
            [imp_midpoint_global_error(i), ~, ~] = implicit_midpoint_global_error(tspan, h(1,i));


            h_avg(1, i) = h_avg_temp;
        end

        %[p1,k1] = loglog_fit(h,analytical_difference);
        [p2,k2] = loglog_fit(h_avg,exp_euler_global_error);
        [p3,k3] = loglog_fit(h_avg,exp_midpoint_global_error);
        [p4,k4] = loglog_fit(h_avg,imp_euler_global_error);
        [p5,k5] = loglog_fit(h_avg,imp_midpoint_global_error);

        % plot the local truncation error
        hold off
        %loglog(h,analytical_difference, 'b.', MarkerSize=10)
        loglog(h_avg,exp_euler_global_error, 'r.', MarkerSize=10)
        hold on
        loglog(h_avg,exp_midpoint_global_error, 'g.', MarkerSize=10)
        loglog(h_avg,imp_euler_global_error, 'c.', MarkerSize=10)
        loglog(h_avg,imp_midpoint_global_error, 'm.', MarkerSize=10)
        %loglog(h, k1*(h.^p1), 'b-')
        %loglog(h_avg, k2*(h_avg.^p2), 'r-')
        %loglog(h_avg, k3*(h_avg.^p3), 'g-')
        %loglog(h_avg, k4*(h_avg.^p4), 'c-')
        %loglog(h_avg, k5*(h_avg.^p5), 'm-')
        xlabel("Average timestep length h")
        ylabel("Error")
        lgd = legend("forward euler global error", "explicit midpoint global error", "backward euler global error", "implicit midpoint global error");
        lgd.Location = "southeast";
        title("Global Truncation Error of All Methods")
    end

    % global error scaling with number of function calls for explicit
    % methods
    if mode == 10
        tspan = [0, pi/6];
        h = logspace(-5, 1, 100);
        exp_euler_global_error = zeros(1, length(h));
        exp_midpoint_global_error = zeros(1, length(h));
        exp_midpoint_evals = zeros(1, length(h));
        exp_euler_evals = zeros(1, length(h));

        for i = 1:length(h)
            [exp_euler_global_error(i), ~, exp_euler_temp] = forward_euler_global_error(tspan, h(1,i));
            [exp_midpoint_global_error(i), ~, exp_midpoint_temp] = explicit_midpoint_global_error(tspan, h(1,i));
            exp_midpoint_evals(i) = exp_midpoint_temp;
            exp_euler_evals(i) = exp_euler_temp;
        end

        [p1,k1] = loglog_fit(exp_euler_evals,exp_euler_global_error);
        [p2,k2] = loglog_fit(exp_midpoint_evals,exp_midpoint_global_error);

        % plot the global truncation error
        hold off
        loglog(exp_euler_evals,exp_euler_global_error, 'r.', MarkerSize=10)
        hold on
        loglog(exp_midpoint_evals,exp_midpoint_global_error, 'g.', MarkerSize=10)
        loglog(exp_euler_evals, k1*(exp_euler_evals.^p1), 'r-')
        loglog(exp_midpoint_evals, k2*(exp_midpoint_evals.^p2), 'g-')
        lgd = legend("forward euler", "explicit midpoint");
        lgd.Location = "southeast";
        xlabel("# of Rate Function Calls")
        ylabel("Error")
        title("Global Truncation Error given # of Function Calls for Explicit Methods")
    end

    % global error compared to number of function calls for all functions
    if mode == 11
        tspan = [0, pi/6];
        h = logspace(-5, 1, 100);
        exp_euler_global_error = zeros(1, length(h));
        exp_midpoint_global_error = zeros(1, length(h));
        imp_euler_global_error = zeros(1, length(h));
        imp_midpoint_global_error = zeros(1, length(h));


        exp_midpoint_evals = zeros(1, length(h));
        exp_euler_evals = zeros(1, length(h));
        imp_midpoint_evals = zeros(1, length(h));
        imp_euler_evals = zeros(1, length(h));


        for i = 1:length(h)
            [exp_euler_global_error(i), ~, exp_euler_temp] = forward_euler_global_error(tspan, h(1,i));
            [exp_midpoint_global_error(i), ~, exp_midpoint_temp] = explicit_midpoint_global_error(tspan, h(1,i));
    
            [imp_euler_global_error(i), ~, imp_euler_temp] = backward_euler_global_error(tspan, h(1,i));
            [imp_midpoint_global_error(i), ~, imp_midpoint_temp] = implicit_midpoint_global_error(tspan, h(1,i));

            exp_midpoint_evals(i) = exp_midpoint_temp;
            exp_euler_evals(i) = exp_euler_temp;
            imp_midpoint_evals(i) = imp_midpoint_temp;
            imp_euler_evals(i) = imp_euler_temp;
        end

        [p1,k1] = loglog_fit(exp_euler_evals,exp_euler_global_error);
        [p2,k2] = loglog_fit(exp_midpoint_evals,exp_midpoint_global_error);
        [p3,k3] = loglog_fit(imp_euler_evals,imp_euler_global_error);
        [p4,k4] = loglog_fit(imp_midpoint_evals,imp_midpoint_global_error);

        % plot the local truncation error
        hold off
        loglog(exp_euler_evals,exp_euler_global_error, 'r.', MarkerSize=10)
        hold on
        loglog(exp_midpoint_evals,exp_midpoint_global_error, 'g.', MarkerSize=10)
        loglog(imp_euler_evals,imp_euler_global_error, 'b.', MarkerSize=10)
        loglog(imp_midpoint_evals,imp_midpoint_global_error, 'm.', MarkerSize=10)
        %loglog(exp_euler_evals, k1*(exp_euler_evals.^p1), 'r-')
        %loglog(exp_midpoint_evals, k2*(exp_midpoint_evals.^p2), 'g-')
        lgd = legend("forward euler", "explicit midpoint", "backward euler", "implicit midpoint");
        lgd.Location = "southeast";
        xlabel("# of Rate Function Calls")
        ylabel("Error")
        title("Global Truncation Error given # of Function Calls for all Methods")
    end

    if mode == 12
        % stability visualization
        % Comparing closed form approximation (of various time step) with numerical solutions
        % Forward Euler Method
        h = 0.45; % WHAT ARE GOOD H VALUES
        tspan = [0, 20];

        % calculate analytical values
        t_analytical_list = linspace(tspan(1), tspan(2), 50);
        analytical_X_list = solution01(t_analytical_list);

        % plot the approximations (of various time step)
        [~, foward_euler_X_list, ~, ~] = forward_euler(@rate_func01,[tspan(1), tspan(2)],solution01(tspan(1)), h);
        [~, explicit_midpoint_X_list, ~, ~] = explicit_midpoint(@rate_func01,[tspan(1), tspan(2)],solution01(tspan(1)), h);
        [~, backward_euler_X_list, ~, ~] = backward_euler(@rate_func01,[tspan(1), tspan(2)],solution01(tspan(1)), h);
        [t_list, implicit_midpoint_X_list, ~, ~] = implicit_midpoint(@rate_func01,[tspan(1), tspan(2)],solution01(tspan(1)), h);
        
        %plot(t_list, foward_euler_X_list, "-")
        
        fig = figure;
        subplot(2,2,1); % 2 row, 2 columns, select the first position
        hold on
        plot(t_analytical_list, analytical_X_list, 'black')
        plot(t_list, foward_euler_X_list, "-")
        title("Forward Euler")

        
        subplot(2,2,2); % 2 row, 2 columns, select the second position
        hold on
        plot(t_analytical_list, analytical_X_list, 'black')
        plot(t_list, explicit_midpoint_X_list, "-")
        title("Explicit Midpoint")

        subplot(2,2,3); % 2 row, 2 columns, select the first position
        hold on
        plot(t_analytical_list, analytical_X_list, 'black')
        plot(t_list, backward_euler_X_list, "-")
        title("Backward Euler")

        subplot(2,2,4); % 2 row, 2 columns, select the second position
        hold on
        plot(t_analytical_list, analytical_X_list, 'black')
        plot(t_list, implicit_midpoint_X_list, "-")
        title("Implicit Midpoint")
        
        han=axes(fig,'visible','off'); 
        han.Title.Visible='on';
        han.XLabel.Visible='on';
        han.YLabel.Visible='on';
        ylabel(han,'X(t)');
        xlabel(han,'t');
        title(han,'Stability Visualizations');
    end
end


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
function dXdt = rate_func01(t,X)
    dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
    X = cos(t);
end

% test func 2
function dXdt = rate_func02(t,X)
    dXdt = [-1*X(2); X(1)];
end

function X = solution02(t)
    X = [cos(t);sin(t)];
end