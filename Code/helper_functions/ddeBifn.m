function bifn = ddeBifn(p)
    % Calculate k1, k2, and k3 based on u* and v*
    kappa.k1 = p.a * p.beta * p.u * (1 - p.u);
    kappa.k2 = p.d * p.beta * p.v * (1 - p.v);
    kappa.k3 = p.b * p.c * p.beta^2 * p.u * p.v * (1 - p.u) * (1 - p.v);
    
    % Preallocate storage for bifurcation points
    bifn.sols = arrayfun(@(x) nan(2, 121), 1:3, 'UniformOutput', false);
    
    % Iterate over tau1 values
    for i = 1:length(p.tau1_vals)
        tau1 = p.tau1_vals(i);
        sols = [];

        % Iterate over initial guesses for tau2 and omega
        for tau2 = p.tau2_vals
            for omega = p.omega_vals
                try
                    % Solve the nonlinear system
                    [x,~,exitflag] = fsolve(@epsilon, [tau2, omega], ...
                        optimoptions('fsolve', 'Display', 'off'));
                    sol = round(x, 3);
    
                    % Verify solution and add if valid
                    if exitflag > 0 && sol(1) >= 0 && sol(1) <= 1.5
                        sols = [sols; sol]; %#ok<AGROW>
                    end
                catch
                    % Ignore failures in solving
                end
            end
        end
        
        % If solutions exist, find unique tau2 values
        if ~isempty(sols)        
            % Remove repeated solutions of tau2
            sols_unique = unique(sols(:,1),'rows');
    
            if tau1 <= 4.6
                if size(sols_unique,1) == 1
                    bifn.sols{3}(:,i) = [tau1, sols_unique(1,1)];
                elseif size(sols_unique,1) == 2
                    bifn.sols{1}(:,i) = [tau1, sols_unique(1,1)];
                    bifn.sols{2}(:,i) = [tau1, sols_unique(2,1)];
                elseif size(sols_unique,1) >= 3
                    bifn.sols{1}(:,i) = [tau1, sols_unique(1,1)];
                    bifn.sols{2}(:,i) = [tau1, sols_unique(2,1)];
                    bifn.sols{3}(:,i) = [tau1, sols_unique(3,1)];
                end
            else
                bifn.sols{1}(:,i) = [tau1, sols_unique(1,1)];
            end
        end
    end
    
    bifn.line_1 = bifn.sols{3};
    bifn.line_1(:, any(isnan(bifn.line_1), 1)) = [];
    
    bifn.line_2 = [fliplr(bifn.sols{2}), bifn.sols{1}];
    bifn.line_2(:, any(isnan(bifn.line_2), 1)) = [];

%% --------------------------------------------------------------------- %%
% ------------------------------- bias(v) ------------------------------- %
    % Define the bias functions for theta_u and theta_v
    function e = epsilon(x)
        e = [
            ReEpsilon(kappa, struct('tau1', tau1, 'tau2', x(1)), p, x(2));
            ImEpsilon(kappa, struct('tau1', tau1, 'tau2', x(1)), p, x(2));
        ];
    end

% ---------------------------- Re ε(iω) = 0 ----------------------------- %
    % Define the pair of equations for θ_u and θ_v
    function ReE = ReEpsilon(kappa, tau, p, omega)
        ReE = (1 - kappa.k1 .* cos(omega .* tau.tau1)) .* ...
            (1 - kappa.k2 .* cos (omega .* tau.tau1)) - ...
            (omega + kappa.k1 .* sin(omega .* tau.tau1)) .* ...
            ((omega / p.alpha) + kappa.k2 .* sin(omega .* tau.tau1)) - ...
            (kappa.k3 .* cos(2 .* omega .* tau.tau2));
    end

% ---------------------------- Im ε(iω) = 0 ----------------------------- %
    % Define the pair of equations for θ_u and θ_v
    function ImE = ImEpsilon(kappa, tau, p, omega)
        ImE = (1 - kappa.k1 .* cos(omega .* tau.tau1)) .* ...
            ((omega / p.alpha) + kappa.k2 .* sin(omega .* tau.tau1)) + ...
            (omega + kappa.k1 .* sin(omega .* tau.tau1)) .* ...
            (1 - kappa.k2 .* cos(omega .* tau.tau1)) + ...
            (kappa.k3 .* sin(2 .* omega .* tau.tau2));
    end
end