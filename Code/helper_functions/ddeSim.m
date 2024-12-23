%% DDE Simulation
function [sol, null] = ddeSim(p)    
    % Calculate nullclines
    null.u_range = linspace(0, 1, 1000);
    null.v_range = linspace(0, 1, 1000);
    null.nullclines = getNullclines(null.u_range, null.v_range, p);
    
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

% ----------------------------- f_inv(x,p) ------------------------------ %
    % Define the inverse sigmoid function, for z = u,v
    function f = f_inv(z,p)
        f = (1 ./ p.beta) .* log(z ./ (1 - z));
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

% ---------------------------- getNullclines(t,y,Z) ---------------------------- %
    function nullclines = getNullclines(u_range, v_range, p)
        nullclines.u = (f_inv(v_range,p) - p.theta_v - (p.d .* v_range)) ./ p.c;
        nullclines.v = (f_inv(u_range,p) - p.theta_u - (p.a .* u_range)) ./ p.b;
    end