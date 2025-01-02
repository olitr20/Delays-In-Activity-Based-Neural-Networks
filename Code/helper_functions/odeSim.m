function [stst_u, stst_v, sol, nullclines] = odeSim(p)
% ODESIM  Simulate a system of ordinary differential equations.
%   Input:
%       p:  structure containing model parameters of the form {\alpha,
%           \beta, a, b, c, d, \theta_{u}, \theta_{v}}.
%   Output:
%       stst_u: u value of the steady state, if present, found by analysing
%           the final portion of xode for minimal movement.
%       stst_v: v value of the steady state, if present, found by analysing
%           the final portion of xode for minimal movement.
%       sol: solutions of the ordinary differential simulation over a
%           specified period of time.
%       nullclines: structure containing calculated u- and v-nullclines
%           along with corresponding ranges for u and v.

    % Define simulation parameters
    x0 = [0.6; 0]; % initial steady state guess
    tspan = [0 60]; % simulation time span
    
    % Simulate system
    [~,sol] = ode23s(@(t, x) odefun(t, x, p), tspan, x0);

    % Extract steady states
    if round(sol(end,:),4) == round(sol(end-1,:),4)
        stst_u = sol(end,1);
        stst_v = sol(end,2);
    else
        fprintf(['odeSim Warning: No steady state found for parameters: \n' ...
            '{alpha=%.0f, beta=%.0f, a=%.0f, b=%.0f, c=%.0f, d=%.0f, theta_u=%.1f, theta_v=%.1f}. \n'], ...
            p.alpha, p.beta, p.a, p.b, p.c, p.d, p.theta_u, p.theta_v)
        stst_u = [];
        stst_v = [];
    end

    % Calculate nullclines
    u_range = linspace(0, 1, 1000);
    v_range = linspace(0, 1, 1000);
    [u_null, v_null] = getNullclines(u_range, v_range, p);
    
    % Save nullclines
    nullclines.u_range = u_range; nullclines.v_range = v_range;
    nullclines.u_null = u_null; nullclines.v_null = v_null;

%% --------------------------------------------------------------------- %%
% ------------------------------- f(z,p) -------------------------------- %
    % Define the inverse sigmoid function, for z = u,v
    function f = f(z,p)
        f = 1 ./ (1 + exp(-p.beta * z));
    end

% ----------------------------- f_inv(z,p) ------------------------------ %
    % Define the inverse sigmoid function, for z = u,v
    function f = f_inv(z,p)
        f = (1 ./ p.beta) .* log(z ./ (1 - z));
    end

% ---------------------------- odefun(t,y,Z) ---------------------------- %
    function d = odefun(~,x,p)
        u = x(1);
        v = x(2);

        dudt = -u + f(p.theta_u + p.a .* u + p.b .* v, p);
        dvdt = p.alpha .* (-v + f(p.theta_v + p.c .* u + p.d .* v, p));
    
        d = [dudt; dvdt];
    end

% ----------------- getNullclines(u_range, v_range, p) ------------------ %
    function [nulls_u, nulls_v] = getNullclines(u_range, v_range, p)
        nulls_u = (f_inv(v_range,p) - p.theta_v - (p.d .* v_range)) ./ p.c;
        nulls_v = (f_inv(u_range,p) - p.theta_u - (p.a .* u_range)) ./ p.b;
    end
end