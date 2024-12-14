%% Figure 1 Replicate
% Select Parameters
p.alpha = 1; p.beta = 1;
p.a = 10; p.b = -10;
p.c = 10; p.d = 2;

% Define u ranges for bifurcation search
hopf.us = [0.18:0.0001:0.2764, 0.7236:0.0001:0.82];
sn.us = 0.1127:0.0001:0.8873;

% Calculate hopf, saddle-node and bogdanov-takens bifurcations
[hopf, sn] = odeBifn(hopf, sn, p);

% Initialise figure 1
figure(1);
clf; hold on;

% Plot hopf bifurcation
plot(hopf.theta.uP, hopf.theta.vP, 'k--', 'LineWidth',1); % plus values
plot(hopf.theta.uM, hopf.theta.vM,'k--', 'LineWidth',1); % minus values

% Plot saddle node bifurcation
plot(sn.theta.uP, sn.theta.vP, 'k', 'LineWidth',1.5); % plus values
plot(sn.theta.uM, sn.theta.vM, 'k', 'LineWidth',1.5); % minus values

% Format axes
xlabel("$\theta_{\mathit{u}}$", 'Interpreter', 'latex')
ylabel("$\theta_{\mathit{v}}$", 'Rotation', 0, 'Interpreter', 'latex')
set(gca,'FontSize', 14, 'FontName', 'Times')
xlim([-7 7]);
ylim([-12 0]);
xticks(-6:2:6)
xticklabels({'','-4','','0','','-4',''})
yticks(-12:2:0)
yticklabels({'-12','','-8','','-4','','0'})

% Graph text
text(0, -9.2, "HB", 'FontSize', 14, 'FontName', 'Times')
text(5, -4, "SN", 'FontSize', 14, 'FontName', 'Times')

% Save figure
print(gcf, '../Figures/Figure_1a.png', '-dpng', '-r300');

%% Figure 1 Replicate - With Bogdanov-Takens Bifurcations
% Select Parameters
p.alpha = 1; p.beta = 1;
p.a = 10; p.b = -10;
p.c = 10; p.d = 2;

% Define u ranges for bifurcation search
hopf.us = [0.18:0.0001:0.2764, 0.7236:0.0001:0.82];
sn.us = 0.1127:0.0001:0.8873;
bt.us = [0.25769, 0.74231];

% Calculate hopf, saddle-node and bogdanov-takens bifurcations
[hopf, sn, bt] = odeBifn(hopf, sn, bt, p);

% Initialise figure 1
figure(1);
clf; hold on;

% Plot hopf bifurcation
plot(hopf.theta.uP, hopf.theta.vP, 'k--', 'LineWidth',1); % plus values
plot(hopf.theta.uM, hopf.theta.vM,'k--', 'LineWidth',1); % minus values

% Plot saddle node bifurcation
plot(sn.theta.uP, sn.theta.vP, 'k', 'LineWidth',1.5); % plus values
plot(sn.theta.uM, sn.theta.vM, 'k', 'LineWidth',1.5); % minus values

% Plot bogdanov-takens bifurcations
plot(bt.theta.uP, bt.theta.vP, 'ro', 'LineWidth',1.5,'MarkerSize',12); % plus values
plot(bt.theta.uM, bt.theta.vM,'ro', 'LineWidth',1.5,'MarkerSize',12); % minus values

% Format axes
xlabel("$\theta_{\mathit{u}}$", 'Interpreter', 'latex')
ylabel("$\theta_{\mathit{v}}$", 'Rotation', 0, 'Interpreter', 'latex')
set(gca,'FontSize', 14, 'FontName', 'Times')
xlim([-7 7]);
ylim([-12 0]);
xticks(-6:2:6)
xticklabels({'','-4','','0','','-4',''})
yticks(-12:2:0)
yticklabels({'-12','','-8','','-4','','0'})

% Graph Text
text(0, -9.2, "HB", 'FontSize', 14, 'FontName', 'Times')
text(5, -4, "SN", 'FontSize', 14, 'FontName', 'Times')

% Save figure
print(gcf, '../Figures/Figure_1b.png', '-dpng', '-r300');

%% Figure 2 Replicate
% Select Parameters
p.alpha = 1; p.beta = 1;
p.a = 10; p.b = -10;
p.c = 10; p.d = 2;

% Find u* and v* such that θu = -2 & θv = -4
p.theta_u = -2;
p.theta_v = -4;
[p.u, p.v] = calcBias(p, p.theta_u, p.theta_v);

% Specify ranges for tau1, tau2, and omega
p.tau1_vals = linspace(0, 12, 121); % Range of tau1
p.tau2_vals = linspace(0, 1.5, 3); % Range of tau2
p.omega_vals = linspace(0, 2.5, 7); % Omega range

% Calculate hopf bifurcation lines in \tau_{1} and \tau_{2}
bifn = ddeBifn(p);

% Initialise figure 2
figure(2);
clf; hold on;

% Plot bifurcation lines
plot(bifn.line_1(1,:), bifn.line_1(2,:), 'k')
plot(bifn.line_2(1,:), bifn.line_2(2,:), 'k')

% Format Axes
xlabel("\tau_1");
ylabel("\tau_2", 'rotation', 0);
set(gca,'FontSize', 14, 'FontName', 'Times')
xlim([0 12]);
ylim([0 1.5]);
xticks([0 4 8 12])
yticks([0 0.5 1.0 1.5])
yticklabels({'0', '0.5', '1.0', '1.5'});

% Add annotations
annotation('textbox',[0.2046 0.5762 0.1219 0.0619],'String','unstable',...
    'Rotation',71.5,'FontSize',14,'FontName','Times','EdgeColor','none');
annotation('textbox',[0.3028 0.5262 0.0969 0.0619],'String','stable',...
    'Rotation',71.5,'FontSize',14,'FontName','Times','EdgeColor','none');
annotation('textbox',[0.5081 0.3905 0.1219 0.0619],'String','unstable',...
    'FontSize',14,'FontName','Times','EdgeColor','none');

% Save figure
print(gcf, '../Figures/Figure_2.png', '-dpng', '-r300');

%% Figure 3 Replicate
% Select Parameters
p.alpha = 1; p.beta = 700; % high beta approximates heaviside
p.a = -1; p.b = -0.4;
p.c = -0.4; p.d = -1;
p.tau1 = 1; p.tau2 = 1.4;

% Find u* and v* such that θu = 0.7 & θv = 0.7
p.theta_u = 0.7;
p.theta_v = 0.7;
[p.u, p.v] = calcBias(p, p.theta_u, p.theta_v);

% Define DDE parameters
p.tspan = [0 15];
p.delays = [p.tau1 p.tau2];
p.history1 = [0.37 0.8];
p.history2 = [0.2 0.8];
p.options = ddeset('RelTol', 1e-5);

% Solve DDE
sol.Synchronous = dde23(@(t, y, Z) ddefun(t, y, Z, p), ...
    p.delays,p.history1,p.tspan,p.options);
sol.AntiSynchronous = dde23(@(t, y, Z) ddefun(t, y, Z, p), ...
    p.delays,p.history2,p.tspan,p.options);

% Initilaise figure 3
figure(3);
clf;
tiledlayout(2,1,'TileSpacing','Compact','Padding','loose');

% First subplot
nexttile
hold on;

% Plot synchronous solution
plot(sol.Synchronous.x, sol.Synchronous.y(1,:), ...
    'k', 'linewidth', 1.5); % u against time
dashline(sol.Synchronous.x,sol.Synchronous.y(2,:), ...
    1.5, 1, 1.5, 1, 'color', '#808080', 'linewidth', 1.5) % v against time

% Format axes
set(gca,'FontSize', 14, 'FontName', 'Times')
% xlim([0, 15]);
ylim([0, 1]);
% xticks([0 5 10 15])
xticklabels({'','','',''})
yticks([0 0.5 1.0])
yticklabels({'0', '0.5', '1.0'});

% Second subplot
nexttile
hold on;

% Plot antisynchronous solution
plot(sol.AntiSynchronous.x, sol.AntiSynchronous.y(1,:), ...
    'k', 'linewidth', 1.5); % u against time
dashline(sol.AntiSynchronous.x,sol.AntiSynchronous.y(2,:), ...
    1.5, 1, 1.5, 1, 'color', '#808080', 'linewidth', 1.5) % v against time

% Format axes
xlabel("$\mathit{t}$", 'Interpreter', 'latex')
set(gca,'FontSize', 14, 'FontName', 'Times')
% xlim([0, 15]);
ylim([0, 1]);
% xticks([0 5 10 15])
yticks([0 0.5 1.0])
yticklabels({'0', '0.5', '1.0'});

% Add annotations
annotation('textbox',[0.15 0.83 0.0446 0.0536],'String','$\mathit{u}$',...
    'FontSize',14,'EdgeColor','none','Interpreter','latex');
annotation('textbox',[0.15 0.63 0.0446 0.0536],'String','$\mathit{v}$',...
    'FontSize',14,'EdgeColor','none','Interpreter','latex');
annotation('textbox',[0.15 0.38 0.0446 0.0536],'String','$\mathit{u}$',...
    'FontSize',14,'EdgeColor','none','Interpreter','latex');
annotation('textbox',[0.15 0.21 0.0446 0.0536],'String','$\mathit{v}$',...
    'FontSize',14,'EdgeColor','none','Interpreter','latex');

% Save figure
print(gcf, '../Figures/Figure_3.png', '-dpng', '-r300');

%% Figure 9a Replicate
% Select Parameters
p.alpha = 1; p.beta = 60;
p.a = -6; p.b = 2.5;
p.c = 2.5; p.d = -6;
p.tau1 = 0.1; p.tau2 = 0.1;

% Find u* and v* such that θu = 0.2 & θv = 0.2
p.theta_u = 0.2;
p.theta_v = 0.2;
[p.u, p.v] = calcBias(p, p.theta_u, p.theta_v);

% Define DDE parameters
p.tspan = [0 40];
p.delays = [p.tau1 p.tau2];
p.history = [0.074 0.077];
p.options = ddeset('RelTol', 1e-5);

% Solve DDE
sol.Chaotic = dde23(@(t, y, Z) ddefun(t, y, Z, p), ...
    p.delays,p.history,p.tspan,p.options);

% Initialise figure 4
figure(8)
clf; hold on

% Plot chotic solution
plot(sol.Chaotic.y(1,:), sol.Chaotic.y(2,:), ...
    'k', 'linewidth', 1.5); % u against v

% Format axes
xlabel("$\mathit{u}$", 'Interpreter', 'latex')
ylabel("$\mathit{v}$", 'Interpreter', 'latex','rotation',0)
set(gca,'FontSize', 14, 'FontName', 'Times')
xlim([0.062, 0.178]);
ylim([0.06, 0.18]);
xticks([0.08 0.1 0.12 0.14 0.16])
yticks([0.06 0.08 0.1 0.12 0.14 0.16 0.18])
xticklabels({'0.08', '', '0.12', '', '0.16', ''});
yticklabels({'0.06', '', '0.10', '', '0.14', '', '0.18'});

% Save figure
print(gcf, '../Figures/Figure_9a.png', '-dpng', '-r300');

%% Figure 9b Replicate
% Select Parameters
p.alpha = 1; p.beta = 40;
p.a = -6; p.b = 2.5;
p.c = 2.5; p.d = -6;
p.tau1 = 0.1; p.tau2 = 0.1;

% Find u* and v* such that θu = 0.2 & θv = 0.2
p.theta_u = 0.2;
p.theta_v = 0.2;
[p.u, p.v] = calcBias(p, p.theta_u, p.theta_v);

% Define DDE parameters
p.tspan = [0 40];
p.delays = [p.tau1 p.tau2];
p.history = [0.074 0.077];
p.options = ddeset('RelTol', 1e-5);

% Solve DDE
sol.QuasiPeriodic = dde23(@(t, y, Z) ddefun(t, y, Z, p), ...
    p.delays,p.history,p.tspan,p.options);

% Initialise figure 4
figure(9)
clf; hold on

% Plot chotic solution
plot(sol.QuasiPeriodic.y(1,:), sol.QuasiPeriodic.y(2,:), ...
    'k', 'linewidth', 1.5); % u against v

% Format axes
xlabel("$\mathit{u}$", 'Interpreter', 'latex')
ylabel("$\mathit{v}$", 'Interpreter', 'latex','rotation',0)
set(gca,'FontSize', 14, 'FontName', 'Times')
xlim([0.07, 0.105]);
ylim([0.07, 0.105]);
xticks([0.07 0.075 0.08 0.085 0.09 0.095 0.1 0.105])
yticks([0.07 0.075 0.08 0.085 0.09 0.095 0.1 0.105])
xticklabels({'', '', '0.08', '', '0.09', '', '0.10', ''});
yticklabels({'', '0.075', '', '0.085', '', '0.095', '', '0.105'});

% Save figure
print(gcf, '../Figures/Figure_9b.png', '-dpng', '-r300');

%% Figure 8 Replicate
% Liapunov exponent attempts to compute the maximum Liapunov exponent of
% the current simulation. The method is pretty simplistic but works on the
% examples I have fed it. Given a solution x(t) at a series of points
% t1, ..., tn I perturb the solution at each time point, integrate the
% equation to the next time point, and compute the logarithm of the rate of
% growth. This is averaged over the whole time series to give an
% approximation. The size of the perturbation is determined by the
% numerical parameter ptbn.

% Define Lyaounov exponent calculation parameters
point = [-6 2.5]; % provide model parameters [a b]
ptbn = 0.001; % distance to perturb points along main trajectory
int = 2; % interval to simulate perturbations for

% Define model parameters
p.alpha = 1; p.beta = 60;
p.a = point(1); p.b = point(2);
p.c = point(2); p.d = point(1);
p.theta_u = 0.2; p.theta_v = 0.2;
p.tau1 = 0.1; p.tau2 = 0.1;

% Define DDE parameters
p.tspan = [0 50];
p.delays = [p.tau1 p.tau2];
p.history = [0.074 0.077];
p.options = ddeset('RelTol', 1e-5);

% Calculate maximal lyapunov exponent
lambda_max = lyapunovExponent(@(t, y, Z) ddefun(t, y, Z, p), p, ptbn, int);
fprintf('Maximal Lyapunov Exponent: %.4f\n', lambda_max)

%% Calculate Maximum Lyapunov Exponents over parameter space
% Define Lyaounov exponent calculation parameterssprintf("$\\mathit{u} \\; \\text{(a = %d)}$", a), 'Interpreter', 'latex'
ptbn = 0.001; % distance to perturb points along main trajectory
int = 15; % interval to simulate perturbations for

% Define grid over a, b parameter space
res = 256; % specify resolution of grid
p.arange = linspace(-10, 0, res);
p.brange = linspace(0, 5, res);
[p.a_grid, p.b_grid] = meshgrid(p.arange, p.brange);

% Define model parameters
p.alpha = 1; p.beta = 60;
p.theta_u = 0.2; p.theta_v = 0.2;
p.tau1 = 0.1; p.tau2 = 0.1;

% Define DDE parameters
p.tspan = [0 60];
p.delays = [p.tau1 p.tau2];
p.history = [0.074 0.077];
p.options = ddeset('RelTol', 1e-5);

% Calculate maximal lyapunov exponents
lambda_max = zeros(size(p.a_grid, 1), size(p.a_grid, 2));
for i = 1:size(p.a_grid, 1)
    for j = 1:size(p.a_grid, 2)
        p.a = p.a_grid(i,j); p.b = p.b_grid(i,j);
        p.c = p.b; p.d = p.a;
        lambda_max(i,j) = lyapunovExponent(@(t,y,Z) ddefun(t,y,Z,p), ...
            p, ptbn, int);
    end
    fprintf('Completed a = %.3f\n', p.a)
end

le_img = lambda_max;
le_img(le_img < 0) = 0;

figure(10);
clf; hold on

imagesc(p.a_grid(1, :), p.b_grid(:, 1), le_img);
colormap(flipud(gray));
colorbar;

xlabel("$\mathit{a}$", 'Interpreter', 'latex')
ylabel("$\mathit{b}$", 'Interpreter', 'latex','rotation',0)
set(gca,'FontSize', 14, 'FontName', 'Times')

astep = (p.arange(2) - p.arange(1)) / 2;
bstep = (p.brange(2) - p.brange(1)) / 2;

xlim([-10-astep, 0+astep]);
ylim([0-bstep, 5+bstep]);
xticks([-10 -8 -6 -4 -2 0])
yticks([0 1 2 3 4 5])
xticklabels({'-10', '-8', '-6', '-4', '-2', '0'});
yticklabels({'0', '1', '2', '3', '4', '5'});

% Save figure
print(gcf, '../Figures/Figure_8.png', '-dpng', '-r300');

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