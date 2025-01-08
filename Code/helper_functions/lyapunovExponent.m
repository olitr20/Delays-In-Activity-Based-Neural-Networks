function lambda_max = lyapunovExponent(p, ptbn, int)
% LYAPUNOVEXPONENT calculate the maximal Lyapunov Exponent of a dynamical
% system for a given set of parameters.
%   Inputs:
%       p:  structure containing model parameters of the form {\alpha,
%           \beta, a, b, c, d, \theta_{u}, \theta_{v}, \tau_{1}, \tau_{2},
%           delays, history, timespan, options}.
%       ptbn: size of the perturbation.
%       int: time interval to simulate perturbations for.
%
%   Output:
%       lambda_max - Approximation of the maximal Lyapunov exponent for a
%       given set of parameters passed with the function handle ddefun.

    % Unpack parameters
    tspan = p.tspan;
    delays = p.delays;
    history = p.history;
    options = p.options;

    % Integrate initial condition over time span
    sol = dde23(@(t,y,Z) ddefun(t,y,Z,p), ...
        delays, history, tspan, options);

    % Extract time points from this solution
    timepoints = tspan(1):int:tspan(2);
    N = length(timepoints);

    % Preallocate storage for logarithm of growth rates
    log_growth_rates = zeros(1, N-1);

    % Loop over each time interval
    for i = 1:N-1
        % Calculate current and next states on main trajectory
        t_current = timepoints(i);
        t_next = timepoints(i+1);
        x_current = deval(sol, t_current);
        x_next = deval(sol, t_next);

        % Generate perturbation vector perpendicular to trajectory
        tang_vec = (deval(sol, timepoints(min(i+1, N))) - x_current) / ...
            (timepoints(min(i+1, N)) - t_current); % compute tangent vec
        perp_vec = [-tang_vec(2), tang_vec(1)]; % rotate 90 degrees
        norm_vec = perp_vec / norm(perp_vec); % normalise vec

        % Randomise sign of perturbation
        rng(1) % for reproducibility
        ptbn_r = ptbn * (randi([0, 1]) * 2 - 1);

        % Calculate perturbation vector
        perturbation = ptbn_r * norm_vec;

        % Integrate perturbed solution over interval
        perturbed_sol = dde23(@(t,y,Z) ddefun(t,y,Z,p), ...
            delays, @perturbed_history, [t_current t_next], options);

        % Extract perturbed state at t_next
        x_next_perturbed = deval(perturbed_sol, t_next);

        % Calculate distance between trajectories at t_next
        distance_perturbed = norm(x_next_perturbed - x_next);
        norm_growth = distance_perturbed / norm(perturbation);

        % Calculate the log rate of growth
        log_growth_rates(i) = log(norm_growth);
    end

    % Average log(growth rate) over the whole time series
    lambda_max = mean(log_growth_rates);

%% --------------------------------------------------------------------- %%
%% --------------------------------------------------------------------- %%
% ------------------------------- f(x,p) -------------------------------- %
    % Define the inverse sigmoid function, for z = u,v
    function f = f(z,p)
        f = 1 ./ (1 + exp(-p.beta * z));
    end

% ---------------------------- ddefun(t,y,Z) ---------------------------- %
    function d = ddefun(~,y,Z,p)
        dudt = -y(1) + ...
            f(p.theta_u + p.a .* Z(1,1) + ...
            p.b .* Z(2,2),p);
        dvdt = p.alpha .* ...
            (-y(2) + f(p.theta_v + p.c .* Z(1,2) + ...
            p.d .* Z(2,1),p));
    
        d = [dudt; dvdt];
    end

% ------------------------- perturbed_history(t) ------------------------ %
    % Define the history function
    function h = perturbed_history(t)
        pt = perturbation';
        if t < 0
            h = p.history' + pt;
        elseif t >= 0 && t <= sol.x(end)
            h = deval(sol, t) + pt;
        end
    end
end
