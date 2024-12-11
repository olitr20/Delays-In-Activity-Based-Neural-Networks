function lambda_max = lyapunovExponent(ddefun, p, ptbn)
% lyapunovExponent computes the maximum Lyapunov exponent of a dynamical system.
% 
% Inputs:
%   ddefun - Function handle defining the system of DDEs (dx/dt = f(t, x, x(t-tau))).
%   p      - Structure containing the parameters: tspan, delays, history, and options.
%   ptbn   - Size of the perturbation (scalar).
%
% Output:
%   lambda_max - Approximation of the maximum Lyapunov exponent.

    % Unpack parameters
    tspan = p.tspan;
    delays = p.delays;
    history = p.history;
    options = p.options;

    % Integrate the initial condition over the time span
    sol = dde23(ddefun, delays, history, tspan, options);

    % Extract time points from the solution
    time_points = sol.x;
    main_trajectory = sol.y;
    N = length(time_points);

    % Preallocate storage for logarithm of growth rates
    log_growth_rates = zeros(1, N-1);

    % Prepare figure for plotting
    figure(10);
    clf; hold on;
    plot(main_trajectory(1, :), main_trajectory(2, :), 'k', 'LineWidth', 1.5); % Main trajectory
    xlabel("$\mathit{u}$", 'Interpreter', 'latex')
    ylabel("$\mathit{v}$", 'Interpreter', 'latex','rotation',0)
    set(gca,'FontSize', 14, 'FontName', 'Times')

    % Initialise Lyapunov exponent accumulator
    sum_log_norms = 0;

    % Loop over each time interval
    for i = 1:N-1
        % Current state
        t_current = time_points(i);
        x_t = deval(sol, t_current);

        % Compute the tangent vector (approximated as the derivative)
        dx_dt = (deval(sol, time_points(min(i+1, N))) - x_t) / (time_points(min(i+1, N)) - t_current);

        % Generate a perturbation vector perpendicular to the trajectory
        perp_vector = [-dx_dt(2), dx_dt(1)]; % Rotate 90 degrees
        ptbn_r = ptbn * (randi([0, 1]) * 2 - 1);
        perturbation = ptbn_r * perp_vector;

        % Integrate the perturbed solution to the next time point
        t_next = time_points(i+1);
        t_max = t_next-t_current;
        perturbed_history = sol.y(:,i)' + perturbation; % Adjust the history for perturbed solution
        perturbed_sol = dde23(ddefun, delays, perturbed_history, [0, t_max], options);

        % Plot the perturbation
        plot(perturbed_sol.y(1,1), perturbed_sol.y(2,1), ...
            '.', 'MarkerSize', 6, 'LineWidth', 1, 'color', '#f77e1b');
        plot(perturbed_sol.y(1,:), perturbed_sol.y(2,:), ...
            '-', 'LineWidth', 1, 'color', '#f77e1b');

        % Extract perturbed state at t_next
        x_t_perturbed = deval(perturbed_sol, t_max);

        % Compute the Euclidean norm of the distance between trajectories
        distance_perturbed = norm(x_t_perturbed - x_t);
        norm_growth = distance_perturbed / norm(perturbation);

        % Compute the logarithm of the rate of growth
        log_growth_rates(i) = log(norm_growth);
    end

    % Average log(growth rate) over the whole time series
    lambda_max = mean(log_growth_rates);

    % Finalise plot
    hold off;
end
