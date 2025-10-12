function [t_list,X_list,h_avg, num_evals] = explicit_midpoint(rate_func_in,tspan,X0,h_ref)
    num_of_steps = ceil((tspan(2) - tspan(1)) / h_ref);
    h_avg = (tspan(2) - tspan(1)) / num_of_steps;
    t_list = linspace(tspan(1), tspan(2), num_of_steps+1)';
    X_list = [];
    X_list(:,1) = X0;
    num_evals = 0;

    for i = 2:length(t_list)
        [XB_temp,num_eval_temp] = explicit_midpoint_step(rate_func_in, t_list(i-1), X_list(:, i-1), h_avg);
        X_list(:, i) = XB_temp;
        num_evals = num_evals + num_eval_temp;
    end
end


%This function computes the value of X at the next time step
%using the explicit midpoint approximation
function [XB,num_evals] = explicit_midpoint_step(rate_func_in,t,XA,h)
    XB_half = XA + (h/2) * rate_func_in(t, XA);
    XB = XA +  h* rate_func_in(t+(h/2), XB_half);
    
    num_evals = 2;
end
