function [u, v] = calcBias(p, theta_u, theta_v)
    
    sol = fsolve(@bias, [0.5, 0.5], ...
        optimoptions('fsolve', 'Display', 'off'));
    
    u = sol(1);
    v = sol(2);

%% --------------------------------------------------------------------- %%
% ----------------------------- f_inv(x,p) ------------------------------ %
    % Define the inverse sigmoid function, for z = u,v
    function f = f_inv(z)
        f = (1 ./ p.beta) .* log(z ./ (1 - z));
    end

% ------------------------------- bias(v) ------------------------------- %
    % Define the bias functions for theta_u and theta_v
    function b = bias(x)
        b = [
            f_inv(x(1)) - (p.a .* x(1)) - (p.b .* x(2)) - theta_u;
            f_inv(x(2)) - (p.c .* x(1)) - (p.d .* x(2)) - theta_v;
        ];
    end
end