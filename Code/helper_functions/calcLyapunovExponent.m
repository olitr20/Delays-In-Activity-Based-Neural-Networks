function lambda_max = calcLyapunovExponent(p)
% CALCLYAPUNOVEXPONENT Calculate maximal Lyapunov Exponents by numerically
% integrationg the variational equation (identical to the underlying system
% at each timepoint but incorporating a perturbed history) in parallel with
% the underlying system. Logarithm rate of growth is calculated at the end
% of a specified time span for perturbation simulations which is then
% averaged over the whole time series.
%   Input:
%       p:  structure containing model parameters of the form {\alpha,
%           \beta, \theta_{u}, \theta_{v}, \tau_{1}, \tau_{2}, delays,
%           history, timespan, options, ptbn (perturbation amount), int
%           (time interval to simulate perturbations for), a_grid (grid of
%           a parameter values at a given resolution), b_grid (grid of b
%           parameter values at a given resolution)}.
%   Output:
%       lambda_max: a numeric array containing all maximal lyapunov
%           components over the specified parameter range set by a_grid and
%           b_grid.

    % Calculate maximal lyapunov exponents
    lambda_max = zeros(size(p.a_grid, 1), size(p.a_grid, 2));
    t_full = tic;
    for i = 1:size(p.a_grid, 1)
        t_iter = tic;
        for j = 1:size(p.a_grid, 2)
            p.a = p.a_grid(i,j); p.b = p.b_grid(i,j);
            p.c = p.b; p.d = p.a;
            lambda_max(i,j) = lyapunovExponent(@(t,y,Z) ddefun(t,y,Z,p), ...
                p, p.ptbn, p.int);
        end
        elapsed_iter = toc(t_iter);
        fprintf('Completed b = %.2f, elapsed time: %.1f seconds\n', p.b, elapsed_iter)
    end
    elapsed_full = toc(t_full);
    fprintf('Total elasped time: %.1f seconds\n', elapsed_full)
end

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