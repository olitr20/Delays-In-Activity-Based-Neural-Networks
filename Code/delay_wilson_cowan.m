%% Modelling Delays in an Activity-Based Wilson-Cowan Neural Network
clearvars; clc;
addpath('Code/helper_functions/');
cd('Code/')

%% Figure 1 Replicate
% Select Parameters
p.alpha = 1; p.beta = 1;
p.a = 10; p.b = -10;
p.c = 10; p.d = 2;

% Calculate hopf, saddle-node and bogdanov-takens bifurcations
[hopf, sn] = odeBifn(p);

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

% Calculate hopf, saddle-node and bogdanov-takens bifurcations
[hopf, sn, bt] = odeBifn(p);

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

% Define a custom colour palette
c_palette = customColourPalette(["f77e1b"; "2c70b3"; "c74440"]);

% Define omega range
omega.min = 0; omega.max = max(cellfun(@max, bifn.omegas));

for i = 1:length(bifn.line_1)-1
    % Define each segment
    x_seg = [bifn.line_1(1,i), bifn.line_1(1,i+1)];
    y_seg = [bifn.line_1(2,i), bifn.line_1(2,i+1)];
    
    % Map omega to colormap
    omega.val = bifn.line_1a(i+1);
    omega.colour_idx = round(((omega.val - omega.min) / (omega.max - omega.min)) * (256 - 1)) + 1;
    omega.colour = c_palette(omega.colour_idx, :);
    
    % Plot the segment with the corresponding colour
    plot(x_seg, y_seg, 'Color', omega.colour, 'LineWidth', 2);
end
clear i x_seg y_seg

for i = 1:length(bifn.line_2)-1
    % Define each segment
    x_seg = [bifn.line_2(1,i), bifn.line_2(1,i+1)];
    y_seg = [bifn.line_2(2,i), bifn.line_2(2,i+1)];
    
    % Map frequency to colormap index without normalisation
    omega.val = bifn.line_2a(i+1);
    omega.colour_idx = round(((omega.val - omega.min) / (omega.max - omega.min)) * (256 - 1)) + 1;
    omega.colour = c_palette(omega.colour_idx, :);
    
    % Plot the segment with the corresponding colour
    plot(x_seg, y_seg, 'Color', omega.colour, 'LineWidth', 2);
end
clear i x_seg y_seg

% Format Axes
xlabel("\tau_1");
ylabel("\tau_2", 'rotation', 0);
set(gca,'FontSize', 14, 'FontName', 'Times')
xlim([0 12]);
ylim([0 1.5]);
xticks([0 4 8 12])
yticks([0 0.5 1.0 1.5])
yticklabels({'0', '0.5', '1.0', '1.5'});

% Add a colour bar for reference
colormap(c_palette);
clim([omega.min omega.max]); % Set the colour axis to match your data range
colorbar;
clear omega c_palette

% Add annotations
annotation('textbox',[0.1746 0.5762 0.1219 0.0619],'String','unstable',...
    'Rotation',71.5,'FontSize',14,'FontName','Times','EdgeColor','none');
annotation('textbox',[0.2728 0.5262 0.0969 0.0619],'String','stable',...
    'Rotation',71.5,'FontSize',14,'FontName','Times','EdgeColor','none');
annotation('textbox',[0.4781 0.3905 0.1219 0.0619],'String','unstable',...
    'FontSize',14,'FontName','Times','EdgeColor','none');

% Save figure
print(gcf, '../Figures/Figure_2.png', '-dpng', '-r300');

%% Figure 3 Replicate
% Select Parameters
p.alpha = 1; p.beta = 700; % high \beta approximates heaviside
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
figure(3); clf;
tiledlayout(2,1,'TileSpacing','Compact','Padding','loose');

% First subplot
nexttile; hold on;

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
nexttile; hold on;

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
annotation('textbox',[0.14 0.81 0.0446 0.0536],'String','$\mathit{u}$',...
    'FontSize',14,'EdgeColor','none','Interpreter','latex');
annotation('textbox',[0.15 0.61 0.0446 0.0536],'String','$\mathit{v}$',...
    'FontSize',14,'EdgeColor','none','Interpreter','latex');
annotation('textbox',[0.15 0.37 0.0446 0.0536],'String','$\mathit{u}$',...
    'FontSize',14,'EdgeColor','none','Interpreter','latex');
annotation('textbox',[0.15 0.18 0.0446 0.0536],'String','$\mathit{v}$',...
    'FontSize',14,'EdgeColor','none','Interpreter','latex');

% Save figure
print(gcf, '../Figures/Figure_3.png', '-dpng', '-r300');

%% Figure 4 - DDE Simulations
% Select Parameters
p.alpha = 1; p.beta = 1;
p.a = 10; p.b = -10;
p.c = 10; p.d = 2;
p.theta_u = -2; p.theta_v = -4;

% Define DDE parameters
p.tspan = [0 50];
p.history = [0.6 0.7];
p.options = ddeset('RelTol', 1e-5);

p.tau1 = 0.5; p.tau2 = 1;
p.delays = [p.tau1 p.tau2];
[sol1, null] = ddeSim(p);

p.tau1 = 3; p.tau2 = 1;
p.delays = [p.tau1 p.tau2];
sol2 = ddeSim(p);

p.tau1 = 6; p.tau2 = 1;
p.delays = [p.tau1 p.tau2];
sol3 = ddeSim(p);

p.tau1 = 10; p.tau2 = 1;
p.delays = [p.tau1 p.tau2];
sol4 = ddeSim(p);

% Initialise figure 4
figure(4); clf;
tiledlayout(2, 2, 'tileSpacing', 'compact', 'padding', 'loose');

% First subplot
nexttile; hold on;

% Plot nullclines
plot(null.nullclines.u, null.v_range, '-', 'color', '#c74440', 'linewidth', 1.5)
plot(null.u_range, null.nullclines.v, '-', 'color', '#2c70b3', 'linewidth', 1.5)

% Plot solution
plot(sol1.y(1,1), sol1.y(2,1), '.', 'color', '#378c47', 'markersize', 12)
plot(sol1.y(1,:), sol1.y(2,:), 'color', '#378c47', 'linewidth', 1)

% Format subplot
ylabel("$\mathit{v}$", 'interpreter', 'latex', 'rotation', 0, ...
    'position', [-0.2, 0.5, -1], 'verticalalignment', 'middle');
set(gca,'fontSize', 14, 'fontName', 'times')
xlim([0 1]); ylim([0 1]);
xticks(0:0.2:1); xticklabels({'0','0.2','0.4','0.6','0.8','1.0'});
yticks(0:0.2:1); yticklabels({'0','0.2','0.4','0.6','0.8','1.0'});

% Second subplot
nexttile; hold on;

% Plot nullclines
plot(null.nullclines.u, null.v_range, '-', 'color', '#c74440', 'linewidth', 1.5)
plot(null.u_range, null.nullclines.v, '-', 'color', '#2c70b3', 'linewidth', 1.5)

% Plot solution
plot(sol2.y(1,1), sol2.y(2,1), '.', 'color', '#378c47', 'markersize', 12)
plot(sol2.y(1,:), sol2.y(2,:), 'color', '#378c47', 'linewidth', 1)

% Format subplot
set(gca,'fontSize', 14, 'fontName', 'times')
xlim([0 1]); ylim([0 1]);
xticks(0:0.2:1); xticklabels({'0','0.2','0.4','0.6','0.8','1.0'});
yticks(0:0.2:1); yticklabels({'0','0.2','0.4','0.6','0.8','1.0'});

% Third subplot
nexttile; hold on;

% Plot nullclines
plot(null.nullclines.u, null.v_range, '-', 'color', '#c74440', 'linewidth', 1.5)
plot(null.u_range, null.nullclines.v, '-', 'color', '#2c70b3', 'linewidth', 1.5)

% Plot solution
plot(sol3.y(1,1), sol3.y(2,1), '.', 'color', '#378c47', 'markersize', 12)
plot(sol3.y(1,:), sol3.y(2,:), 'color', '#378c47', 'linewidth', 1)

% Format subplot
xlabel("$\mathit{u}$", 'interpreter', 'latex')
ylabel("$\mathit{v}$", 'interpreter', 'latex', 'rotation', 0, ...
    'position', [-0.2, 0.5, -1], 'verticalalignment', 'middle')
set(gca,'FontSize', 14, 'FontName', 'Times')
xlim([0 1]); ylim([0 1]);
xticks(0:0.2:1); xticklabels({'0','0.2','0.4','0.6','0.8','1.0'});
yticks(0:0.2:1); yticklabels({'0','0.2','0.4','0.6','0.8','1.0'});

% Fourth subplot
nexttile; hold on;

% Plot nullclines
plot(null.nullclines.u, null.v_range, '-', 'color', '#c74440', 'linewidth', 1.5)
plot(null.u_range, null.nullclines.v, '-', 'color', '#2c70b3', 'linewidth', 1.5)

% Plot solution
plot(sol4.y(1,1), sol4.y(2,1), '.', 'color', '#378c47', 'markersize', 12)
plot(sol4.y(1,:), sol4.y(2,:), 'color', '#378c47', 'linewidth', 1)

% Format subplot
xlabel("$\mathit{u}$", 'Interpreter', 'latex')
set(gca,'FontSize', 14, 'FontName', 'Times')
xlim([0 1]); ylim([0 1]);
xticks(0:0.2:1); xticklabels({'0','0.2','0.4','0.6','0.8','1.0'});
yticks(0:0.2:1); yticklabels({'0','0.2','0.4','0.6','0.8','1.0'});

% Save figure
print(gcf, '../Figures/Figure_4.png', '-dpng', '-r300');

%% Figure 6 Replicate
% Define model parameters
p.alpha = 1; p.beta = 60;
p.a = -1; p.b = -0.4;
p.c = -1; p.d = 0;
p.theta_u = 0.7; p.theta_v = 0.5;
p.tau_1 = 0.2; p.tau_2 = p.tau_1;

bifn = ddeBiftoolMain(p);

% Initialise figure 6
figure(6); clf;
hold on

% Plot bifurcations
plot(bifn.sn1.x, bifn.sn1.y, 'k', 'linewidth', 1) % first saddle-node
plot(bifn.sn2.x, bifn.sn2.y, 'k', 'linewidth', 1) % second saddle-node
plot(bifn.hopf.x, bifn.hopf.y, 'k--', 'linewidth', 1) % hopf
plot(bifn.snpo1.x, bifn.snpo1.y, 'k-o', 'markersize', 4, 'linewidth', 1) % first saddle-node of periodic orbits
plot(bifn.snpo2.x, bifn.snpo2.y, 'k-o', 'markersize', 4, 'linewidth', 1) % second saddle-node of periodic orbits

% Format figure 6
set(gca,'FontSize', 14, 'FontName', 'Times')
xlabel("$\theta_{\mathit{u}}$", 'interpreter', 'latex')
ylabel('\tau', 'rotation', 0, ...
    'position', [0.34, 0.25, -1], 'verticalalignment', 'middle');
xlim([0.4, 1]);
ylim([0, 0.5]);
xticks([0.4 0.5 0.6 0.7 0.8 0.9 1])
xticklabels({'0.4','0.5','0.6','0.7', '0.8', '0.9', '1.0'})
yticks([0 0.1 0.2 0.3 0.4 0.5])
yticklabels({'0', '0.1', '0.2', '0.3', '0.5', '0.5'});

% Save figure
print(gcf, '../Figures/Figure_6.png', '-dpng', '-r300');

%% Figure 7 Recplicate
% Define model parameters
p.alpha = 1; p.beta = 60;
p.a = -1; p.b = -0.4;
p.c = -1; p.d = 0;
p.theta_u = 0.7; p.theta_v = 0.5;

p.tau_1 = 0.5; p.tau_2 = p.tau_1;
[stst.a, po.a] = ddeBiftoolSection(p);

p.tau_1 = 0.2; p.tau_2 = p.tau_1;
[stst.b, po.b] = ddeBiftoolSection(p);

p.tau_1 = 0.09; p.tau_2 = p.tau_1;
[stst.c, po.c] = ddeBiftoolSection(p);

% Initialise figure 7
figure(7); clf
tiledlayout(3,1,'TileSpacing','Compact','Padding','loose');

% Subplot 1 - \tau = 0.5
nexttile; hold on;

plot(stst.a.stable1.x, stst.a.stable1.y, 'k-')
plot(stst.a.unstable.x, stst.a.unstable.y, 'k--')
plot(stst.a.stable2.x, stst.a.stable2.y, 'k-')

if isfield(po.a, 'stable1')
    plot(po.a.stable1.x, po.a.stable1.y, 'k-o', 'markersize', 4);
end
if isfield(po.a, 'unstable')
    plot(po.a.unstable.x, po.a.unstable.y, 'k-x', 'markersize', 4);
end
if isfield(po.a, 'stable2')
    plot(po.a.stable2.x, po.a.stable2.y, 'k-o', 'markersize', 4);
end

% Format tile 1
ylabel("$\mathit{u}$", 'rotation', 0, 'interpreter', 'latex', ...
    'position', [0.35, 0.5, -1], 'verticalalignment', 'middle');
xlim([0.4 1]);
ylim([0 1]);
xticks([0.4 0.5 0.6 0.7 0.8 0.9 1])
xticklabels({'0.4','0.5','0.6','0.7', '0.8', '0.9', '1.0'})
yticks([0 0.5 1])
yticklabels({'0', '0.5', '1.0'});

% Subplot 2 - \tau = 0.2
nexttile; hold on;

plot(stst.b.stable1.x, stst.b.stable1.y, 'k-')
plot(stst.b.unstable.x, stst.b.unstable.y, 'k--')
plot(stst.b.stable2.x, stst.b.stable2.y, 'k-')

if isfield(po.b, 'stable1')
    plot(po.b.stable1.x, po.b.stable1.y, 'k-o', 'markersize', 4);
end
if isfield(po.b, 'unstable')
    plot(po.b.unstable.x, po.b.unstable.y, 'k-x', 'markersize', 4);
end
if isfield(po.b, 'stable2')
    plot(po.b.stable2.x, po.b.stable2.y, 'k-o', 'markersize', 4);
end

% Format tile 2
ylabel("$\mathit{u}$", 'rotation', 0, 'interpreter', 'latex', ...
    'position', [0.465, 0.45, -1], 'verticalalignment', 'middle');
xlim([0.5 0.9]);
ylim([0 0.9]);
xticks([0.5 0.6 0.7 0.8 0.9])
xticklabels({'0.5','0.6','0.7','0.8', '0.9'})
yticks([0.2 0.4 0.6 0.8])
yticklabels({'0.2', '0.4', '0.6', '0.8'});

% Subplot 3 - \tau = 0.09
nexttile; hold on;

plot(stst.c.stable1.x, stst.c.stable1.y, 'k-')
plot(stst.c.unstable.x, stst.c.unstable.y, 'k--')
plot(stst.c.stable2.x, stst.c.stable2.y, 'k-')

if isfield(po.c, 'stable1')
    plot(po.c.stable1.x, po.c.stable1.y, 'k-o', 'markersize', 4);
end
if isfield(po.c, 'unstable')
    plot(po.c.unstable.x, po.c.unstable.y, 'k-x', 'markersize', 4);
end
if isfield(po.c, 'stable2')
    plot(po.c.stable2.x, po.c.stable2.y, 'k-o', 'markersize', 4);
end

% Format tile 3
xlabel("$\theta_{\mathit{u}}$", 'Interpreter', 'latex');
ylabel("$\mathit{u}$", 'Rotation', 0, 'Interpreter', 'latex', ...
    'position', [0.465, 0.5, -1], 'verticalalignment', 'middle');
xlim([0.5 0.9]);
ylim([0.3 0.7]);
xticks([0.5 0.6 0.7 0.8 0.9])
xticklabels({'0.5','0.6','0.7', '0.8', '0.9'})
yticks([0.3 0.4 0.5 0.6 0.7])
yticklabels({'0.3', '0.4', '0.5', '0.6', '0.7'});

% Save figure
print(gcf, '../Figures/Figure_7.png', '-dpng', '-r300');

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

%% Calculate Maximum Lyapunov Exponents over parameter space
% This section takes approx 2.5 hours to complete

% Define Lyaounov exponent calculation parameters
p.ptbn = 0.001; % distance to perturb points along main trajectory
p.int = 5; % interval to simulate perturbations for

% Define grid over a, b parameter space
p.res = 512; % specify resolution of grid
p.arange = linspace(-10, 0, p.res); % linear spacing over first parameter
p.brange = linspace(0, 5, p.res); % linear spacing over second parameter
[p.a_grid, p.b_grid] = meshgrid(p.arange, p.brange); % mesh grid over parameter space

% Define model parameters
p.alpha = 1; p.beta = 60;
p.theta_u = 0.2; p.theta_v = 0.2;
p.tau1 = 0.1; p.tau2 = 0.1;

% Define DDE parameters
p.tspan = [0 60];
p.delays = [p.tau1 p.tau2];
p.history = [0.074 0.077];
p.options = ddeset('RelTol', 1e-5);

% Calculate maximum Lyapunov Exponents
% lambda_max = calcLyapunovExponent(p);

% Alternatively, load pre-calculated lambda_max
load('lambda_max.mat')
%% Figure 8 Replictae
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

clear astep bstep

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