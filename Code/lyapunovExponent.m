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

    % Initialise figure 10
    figure(10);
    clf; hold on;

    % Plot main trajectory
    plot(main_trajectory(1, :), main_trajectory(2, :), ...
        'k', 'LineWidth', 1.5);

    % Format axes
    xlabel("$\mathit{u}$", 'Interpreter', 'latex')
    ylabel("$\mathit{v}$", 'Interpreter', 'latex','rotation',0)
    set(gca,'FontSize', 14, 'FontName', 'Times')

    % Loop over each time step
    for i = 1:N-1
        % Calculate current and next states on main trajectory
        t_current = time_points(i);
        t_next = time_points(i+1);
        x_current = deval(sol, t_current);
        x_next = deval(sol, t_next);

        % Generate a perturbation vector perpendicular to trajectory
        tang_vec = (deval(sol, time_points(min(i+1, N))) - x_current) / ...
            (time_points(min(i+1, N)) - t_current); % compute tangent vec
        perp_vec = [-tang_vec(2), tang_vec(1)]; % rotate 90 degrees
        norm_vec = perp_vec / norm(perp_vec); % normalise vec

        % Randomise sign of perturbation
        ptbn_r = ptbn * (randi([0, 1]) * 2 - 1);

        % Calculate perturbation vector
        perturbation = ptbn_r * norm_vec;

        % Integrate the perturbed solution to the next time point
        perturbed_sol = dde23(ddefun, delays, @perturbed_history, ...
            [t_current t_next], options);
        % history_pt = deval(sol, t_current) + perturbation';
        % perturbed_sol = dde23(ddefun, delays, history_pt, ...
        %     [t_current t_next], options);

        % Plot the perturbation
        plot(perturbed_sol.y(1,1), perturbed_sol.y(2,1), ...
            '.', 'MarkerSize', 6, 'LineWidth', 1, 'color', '#f77e1b');
        plot(perturbed_sol.y(1,:), perturbed_sol.y(2,:), ...
            '-', 'LineWidth', 1, 'color', '#f77e1b');

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

    % Finalise plot
    hold off;

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

        % if t < t_current
        %     if t < 0
        %         h = p.history'; % original history
        %     else
        %         h = deval(sol, t); % original trajectory
        %     end
        % elseif t == t_current
        %     h = deval(sol, t) + pt; % apply perturbation at t_current
        % else
        %     h = deval(sol, t); % let system evolve naturally beyond t_current
        % end
    end
end
