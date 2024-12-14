function lambda_max = lyapunovExponent(ddefun, p, ptbn, int)
% lyapunovExponent computes the maximum Lyapunov exponent of a dynamical system.
% 
% Inputs:
%   ddefun - Function handle defining the system of DDEs (dx/dt = f(t, x, x(t-tau))).
%   p      - Structure containing the parameters: tspan, delays, history, and options.
%   ptbn   - Size of the perturbation (scalar).
%   int    - Time interval to simulate perturbations for.
%
% Output:
%   lambda_max - Approximation of the maximum Lyapunov exponent.

    % Unpack parameters
    tspan = p.tspan;
    delays = p.delays;
    history = p.history;
    options = p.options;

    % Integrate initial condition over time span
    sol = dde23(ddefun, delays, history, tspan, options);

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
        perturbed_sol = dde23(ddefun, delays, @perturbed_history, ...
            [t_current t_next], options);

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
