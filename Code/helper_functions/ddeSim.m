function [sol, nullclines] = ddeSim(p)
% DDESIM  Simulate a system of delayed differential equations.
%   Input:
%       p:  structure containing model parameters of the form {\alpha,
%           \beta, a, b, c, d, \theta_{u}, \theta_{v}, delays, history,
%           tspan, options}.
%   Output:
%       sol: solutions of the delayed differential simulation over a
%           specified period of time.
%       nullclines: structure containing calculated u- and v-nullclines
%           along with corresponding ranges for u and v.


    % Calculate nullclines
    u_range = linspace(0, 1, 1000);
    v_range = linspace(0, 1, 1000);
    [u_null, v_null] = getNullclines(u_range, v_range, p);
    
    % Save nullclines
    nullclines.u_range = u_range; nullclines.v_range = v_range;
    nullclines.u_null = u_null; nullclines.v_null = v_null;

    % Solve DDE
    sol = dde23(@(t, y, Z) ddefun(t, y, Z, p), ...
        p.delays,p.history,p.tspan,p.options);
end

%% --------------------------------------------------------------------- %%
% ------------------------------- f(x,p) -------------------------------- %
    % Define the inverse sigmoid function, for z = u,v
    function f = f(z,p)
        f = 1 ./ (1 + exp(-p.beta * z));
    end

% ----------------------------- f_inv(z,p) ------------------------------ %
    % Define the inverse sigmoid function, for z = u,v
    function f = f_inv(z,p)
        f = (1 ./ p.beta) .* log(z ./ (1 - z));
    end

% --------------------------- ddefun(~,y,Z,p) --------------------------- %
    function d = ddefun(~,y,Z,p)
        dudt = -y(1) + ...
            f(p.theta_u + p.a .* Z(1,1) + ...
            p.b .* Z(2,2),p);
        dvdt = p.alpha .* ...
            (-y(2) + f(p.theta_v + p.c .* Z(1,2) + ...
            p.d .* Z(2,1),p));
    
        d = [dudt; dvdt];
    end

% ----------------- getNullclines(u_range, v_range, p) ------------------ %
    function [nulls_u, nulls_v] = getNullclines(u_range, v_range, p)
        nulls_u = (f_inv(v_range,p) - p.theta_v - (p.d .* v_range)) ./ p.c;
        nulls_v = (f_inv(u_range,p) - p.theta_u - (p.a .* u_range)) ./ p.b;
    end