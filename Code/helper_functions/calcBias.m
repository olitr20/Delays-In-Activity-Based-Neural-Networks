function [u, v] = calcBias(p)
% CALCBIAS  Calculate u and v values for a given \theta_{u} and \theta_{v}.
%   Input:
%       p:  structure containing model parameters of the form {\alpha,
%           \beta, a, b, c, d, \theta_{u}, \theta_{v}}.
%   Output:
%       u: calculated u value for a specific \theta_{u} and \theta_{v}.
%       v: calculated v value for a specific \theta_{u} and \theta_{v}.
    
    sol = fsolve(@bias, [0.5, 0.5], ...
        optimoptions('fsolve', 'Display', 'off'));
    
    u = sol(1);
    v = sol(2);

%% --------------------------------------------------------------------- %%
% ----------------------------- f_inv(z,p) ------------------------------ %
    % Define the inverse sigmoid function, for z = u,v
    function f = f_inv(z)
        f = (1 ./ p.beta) .* log(z ./ (1 - z));
    end

% ------------------------------- bias(x) ------------------------------- %
    % Define the bias functions for theta_u and theta_v
    function b = bias(x)
        b = [
            f_inv(x(1)) - (p.a .* x(1)) - (p.b .* x(2)) - p.theta_u;
            f_inv(x(2)) - (p.c .* x(1)) - (p.d .* x(2)) - p.theta_v;
        ];
    end
end