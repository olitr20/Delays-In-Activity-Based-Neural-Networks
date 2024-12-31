function [stst_u, stst_v] = odeSim(p)
% ODESIM  Simulate a system of ordinary differential equations.
%   Input:
%       p:  structure containing model parameters of the form {\alpha,
%           \beta, a, b, c, d, \theta_{u}, \theta_{v}}.
%   Output:
%       sol: solutions of the delayed differential simulation over a
%           specified period of time.
%       nullclines: structure containing calculated u- and v-nullclines
%           along with corresponding ranges for u and v.

    % Define simulation parameters
    x0 = [1; 0]; % initial steady state guess
    tspan = [0 60]; % simulation time span
    
    [~,xode] = ode23s(@(t, x) odefun(t, x, p), tspan, x0);
    
    stst_u = xode(end,1);
    stst_v = xode(end,2);

%% --------------------------------------------------------------------- %%
% ------------------------------- f(z,p) -------------------------------- %
    % Define the inverse sigmoid function, for z = u,v
    function f = f(z,p)
        f = 1 ./ (1 + exp(-p.beta * z));
    end

% ---------------------------- odefun(t,y,Z) ---------------------------- %
    function d = odefun(~,x,p)
        u = x(1);
        v = x(2);

        dudt = -u + f(p.theta_u + p.a .* u + p.b .* v, p);
        dvdt = p.alpha .* (-v + f(p.theta_v + p.c .* u + p.d .* v, p));
    
        d = [dudt; dvdt];
    end
end