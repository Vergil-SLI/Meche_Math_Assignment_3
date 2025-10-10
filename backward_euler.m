function [t_list,X_list,h_avg, num_evals] = backward_euler(rate_func_in,tspan,X0,h_ref)
    num_of_steps = ceil((tspan(2) - tspan(1)) / h_ref);
    h_avg = (tspan(2) - tspan(1)) / num_of_steps;
    t_list = linspace(tspan(1), tspan(2), num_of_steps+1)';
    X_list = zeros(length(X0),length(t_list));
    X_list(:,1) = X0;
    num_evals = 0;

    for i = 2:length(t_list)
        [XB_temp,num_eval_temp] = backward_euler_step(rate_func_in, t_list(i-1), X_list(i-1), h_avg);
        X_list(:,i) = XB_temp;
        num_evals = num_evals + num_eval_temp;
    end

end

%This function computes the value of X at the next time step
%using the Backward Euler approximation
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%t: the value of time at the current step
%XA: the value of X(t)
%h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
%OUTPUTS:
%XB: the approximate value for X(t+h) (the next step)
% formula depends on the integration method used
%num_evals: A count of the number of times that you called
% rate_func_in when computing the next step
function [XB,num_evals] = backward_euler_step(rate_func_in,t,XA,h)
    G = @(Xnplus1) XA + h * rate_func_in(t+h, Xnplus1) - Xnplus1;

    params = struct();
    [XB, ~, num_evals] = multi_newton(G, XA, params);
end