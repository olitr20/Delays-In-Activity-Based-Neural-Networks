%% DDE Simulation
% Select Parameters
p.alpha = 1; p.beta = 1;
p.a = 10; p.b = -10;
p.c = 10; p.d = 2;
p.tau1 = 0.5; p.tau2 = 1;
% p.tau1 = 3; p.tau2 = 1;
% p.tau1 = 6; p.tau2 = 1;
p.theta_u = -2; p.theta_v = -4;

% Define DDE parameters
p.tspan = [0 50];
p.delays = [p.tau1 p.tau2];
p.history = [0.6 0.7];
p.options = ddeset('RelTol', 1e-5);

% Calculate nullclines
u_range = linspace(0, 1, 1000);
v_range = linspace(0, 1, 1000);
nullclines = getNullclines(u_range, v_range, p);

% Solve DDE
sol = dde23(@(t, y, Z) ddefun(t, y, Z, p), ...
    p.delays,p.history,p.tspan,p.options);

% Initialise figure 1
figure(1);
clf; hold on

% Plot nullclines
plot(nullclines.u, v_range, '-', 'color', '#c74440', 'linewidth', 1.5)
plot(u_range, nullclines.v, '-', 'color', '#2c70b3', 'linewidth', 1.5)

% Plot solution
plot(sol.y(1,1), sol.y(2,1), '.', 'color', '#378c47', 'markersize', 12)
for i = 1:length(sol.y)-1
    plot(sol.y(1,i:i+1), sol.y(2,i:i+1), ...
        'color', '#378c47', 'linewidth', 1); % u against v
    % drawnow;
end
clear i

% Format axes
xlabel("$\mathit{u}$", 'Interpreter', 'latex')
ylabel("$\mathit{v}$", 'Interpreter', 'latex','rotation',0)
set(gca,'FontSize', 14, 'FontName', 'Times')
xlim([0 1]);
ylim([0 1]);
xticks(0:0.2:1)
xticklabels({'0','0.2','0.4','0.6','0.8','1.0'})
yticks(0:0.2:1)
yticklabels({'0','0.2','0.4','0.6','0.8','1.0'})

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