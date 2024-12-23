%% Calculate Maximum Lyapunov Exponents over parameter space
% This takes approx 2.5 hours to complete

function lambda_max = calcLyapunovExponent(p)    
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