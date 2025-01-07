%% Modelling Delays in an Activity-Based Wilson-Cowan Neural Network
clearvars; clc;
addpath('helper_functions/');

%% Figure 1 Replicate
% Define Model Parameters
p.alpha = 1; p.beta = 1;
p.a = 10; p.b = -10;
p.c = 10; p.d = 2;

% Calculate hopf, saddle-node and bogdanov-takens bifurcations
[hopf, sn] = odeBifn(p);

% Initialise figure 1
figure(1);
clf; hold on;

% Plot hopf bifurcation
plot(hopf.theta.uP, hopf.theta.vP, 'k--', 'linewidth',1); % plus values
plot(hopf.theta.uM, hopf.theta.vM,'k--', 'linewidth',1); % minus values

% Plot saddle node bifurcation
plot(sn.theta.uP, sn.theta.vP, 'k', 'linewidth',1.5); % plus values
plot(sn.theta.uM, sn.theta.vM, 'k', 'linewidth',1.5); % minus values

% Format axes
xlabel("$\theta_{\mathit{u}}$", 'interpreter', 'latex');
ylabel("$\theta_{\mathit{v}}$", 'rotation', 0, 'interpreter', 'latex');

xlim([-7 7]);
xticks(-6:2:6);
xticklabels({'','-4','','0','','-4',''});

ylim([-12 0]);
yticks(-12:2:0);
yticklabels({'-12','','-8','','-4','','0'});

set(gca,'fontsize', 14, 'fontname', 'times');

% Graph text
text(0, -9.2, "HB", 'fontsize', 14, 'fontname', 'times');
text(5, -4, "SN", 'fontsize', 14, 'fontname', 'times');

% Save figure
print(gcf, '../Figures/Figure_1a.png', '-dpng', '-r300');

%% Figure 1 Replicate - With Bogdanov-Takens Bifurcations
% Define Model Parameters
p.alpha = 1; p.beta = 1;
p.a = 10; p.b = -10;
p.c = 10; p.d = 2;

% Calculate hopf, saddle-node and bogdanov-takens bifurcations
[hopf, sn, bt] = odeBifn(p);

% Initialise figure 1
figure(1);
clf; hold on;

% Plot hopf bifurcation
plot(hopf.theta.uP, hopf.theta.vP, 'k--', 'linewidth',1); % plus values
plot(hopf.theta.uM, hopf.theta.vM,'k--', 'linewidth',1); % minus values

% Plot saddle node bifurcation
plot(sn.theta.uP, sn.theta.vP, 'k', 'linewidth',1.5); % plus values
plot(sn.theta.uM, sn.theta.vM, 'k', 'linewidth',1.5); % minus values

% Plot bogdanov-takens bifurcations
plot(bt.theta.uP, bt.theta.vP, 'ro', 'linewidth',1.5,'markersize',12); % plus values
plot(bt.theta.uM, bt.theta.vM,'ro', 'linewidth',1.5,'markersize',12); % minus values

% Plot simulation points
plot(-4.1, -7, '*', 'color', '#c74440', 'markersize', 5);
plot(-3.7, -7, '*', 'color', '#2c70b3', 'markersize', 5);
plot(-3.3, -7, '*', 'color', '#378c47', 'markersize', 5);
plot(-2.9, -7, '*', 'color', '#f77e1b', 'markersize', 5);

% Plot figure 3 point
plot(-2, -4, '*', 'color', 'k', 'markersize', 5);

% Format axes
xlabel("$\theta_{\mathit{u}}$", 'interpreter', 'latex');
ylabel("$\theta_{\mathit{v}}$", 'rotation', 0, 'interpreter', 'latex');

xlim([-7 7]);
xticks(-6:2:6);
xticklabels({'','-4','','0','','-4',''});

ylim([-12 0]);
yticks(-12:2:0);
yticklabels({'-12','','-8','','-4','','0'});

set(gca,'fontsize', 14, 'fontname', 'times');

% Graph Text
text(0, -9.2, "HB", 'fontsize', 14, 'fontname', 'times');
text(5, -4, "SN", 'fontsize', 14, 'fontname', 'times');

% Save figure
print(gcf, '../Figures/Figure_1b.png', '-dpng', '-r300');

%% Simulate ODE System Around Bifurcation Lines
% Define Model Parameters
p.alpha = 1; p.beta = 1;
p.a = 10; p.b = -10;
p.c = 10; p.d = 2;

% Define simulation parameters
p.tspan = [0 60]; p.history = [0.7; 0.1];

% Initialise figure 2
figure(2); clf;
tiledlayout(2,2,'tilespacing','compact','padding','loose');

% Subplot 1
nexttile; hold on;

% Calculate nullclines and simulate
p.theta_u = -4.1; p.theta_v = -7;
[~, ~, sol.ode_sim_1, nullclines] = odeSim(p);

% Plot label
plot(-0.04, 0.73, '*', 'color', '#c74440', 'markersize', 12);
plot(-0.04, 0.73, 'square', 'color', '#c74440', 'markersize', 20);

% Plot nullclines
plot(nullclines.u_null, nullclines.v_range, ...
    '-', 'color', '#c74440', 'linewidth', 1.5);
plot(nullclines.u_range, nullclines.v_null, ...
    '-', 'color', '#2c70b3', 'linewidth', 1.5);

% Plot simulation
plot(sol.ode_sim_1(1,1), sol.ode_sim_1(1,2), ...
    '.', 'color', '#378c47', 'markersize', 8);
plot(sol.ode_sim_1(:,1), sol.ode_sim_1(:,2), ...
    '-', 'color', '#378c47', 'linewidth', 1);
plot(sol.ode_sim_1(end,1), sol.ode_sim_1(end,2), ...
    '.', 'color', 'k', 'markersize', 12);

% Format axes
ylabel("$\mathit{v}$", 'rotation', 0, 'interpreter', 'latex');

xlim([-0.1 0.8]);
ylim([-0.1 0.8]);

set(gca,'fontsize', 14, 'fontname', 'times');

% Subplot 2
nexttile; hold on;

% Calculate nullclines and simulate
p.theta_u = -3.7; p.theta_v = -7;
[~, ~, sol.ode_sim_2, nullclines] = odeSim(p);

% Plot label
plot(-0.04, 0.73, '*', 'color', '#2c70b3', 'markersize', 12);
plot(-0.04, 0.73, 'square', 'color', '#2c70b3', 'markersize', 20);

% Plot nullclines
plot(nullclines.u_null, nullclines.v_range, ...
    '-', 'color', '#c74440', 'linewidth', 1.5);
plot(nullclines.u_range, nullclines.v_null, ...
    '-', 'color', '#2c70b3', 'linewidth', 1.5);

% Plot simulation
plot(sol.ode_sim_2(1,1), sol.ode_sim_2(1,2), ...
    '.', 'color', '#378c47', 'markersize', 8);
plot(sol.ode_sim_2(:,1), sol.ode_sim_2(:,2), ...
    '-', 'color', '#378c47', 'linewidth', 1);
plot(sol.ode_sim_2(end,1), sol.ode_sim_2(end,2), ...
    '.', 'color', 'k', 'markersize', 12);

% Format axes
xlim([-0.1 0.8]);
ylim([-0.1 0.8]);

set(gca,'fontsize', 14, 'fontname', 'times');

% Subplot 3
nexttile; hold on;

% Calculate nullclines and simulate
p.theta_u = -3.3; p.theta_v = -7;
[~, ~, sol.ode_sim_3, nullclines] = odeSim(p);

% Plot label
plot(-0.04, 0.73, '*', 'color', '#378c47', 'markersize', 12);
plot(-0.04, 0.73, 'square', 'color', '#378c47', 'markersize', 20);

% Plot nullclines
plot(nullclines.u_null, nullclines.v_range, ...
    '-', 'color', '#c74440', 'linewidth', 1.5);
plot(nullclines.u_range, nullclines.v_null, ...
    '-', 'color', '#2c70b3', 'linewidth', 1.5);

% Plot simulation
plot(sol.ode_sim_3(1,1), sol.ode_sim_3(1,2), ...
    '.', 'color', '#378c47', 'markersize', 8);
plot(sol.ode_sim_3(:,1), sol.ode_sim_3(:,2), ...
    '-', 'color', '#378c47', 'linewidth', 1);
plot(sol.ode_sim_3(end,1), sol.ode_sim_3(end,2), ...
    '.', 'color', 'k', 'markersize', 12);

% Format axes
xlabel("$\mathit{u}$", 'interpreter', 'latex');
ylabel("$\mathit{v}$", 'rotation', 0, 'interpreter', 'latex');

xlim([-0.1 0.8]);
ylim([-0.1 0.8]);

set(gca,'fontsize', 14, 'fontname', 'times');

% Subplot 4
nexttile; hold on;

% Calculate nullclines and simulate
p.theta_u = -2.9; p.theta_v = -7;
[~, ~, sol.ode_sim_4, nullclines] = odeSim(p);

% Plot label
plot(-0.04, 0.73, '*', 'color', '#f77e1b', 'markersize', 12);
plot(-0.04, 0.73, 'square', 'color', '#f77e1b', 'markersize', 20);

% Plot nullclines
plot(nullclines.u_null, nullclines.v_range, ...
    '-', 'color', '#c74440', 'linewidth', 1.5);
plot(nullclines.u_range, nullclines.v_null, ...
    '-', 'color', '#2c70b3', 'linewidth', 1.5);

% Plot simulation
plot(sol.ode_sim_4(1,1), sol.ode_sim_4(1,2), ...
    '.', 'color', '#378c47', 'markersize', 8);
plot(sol.ode_sim_4(:,1), sol.ode_sim_4(:,2), ...
    '-', 'color', '#378c47', 'linewidth', 1);

% Format axes
xlabel("$\mathit{u}$", 'interpreter', 'latex');

xlim([-0.1 0.8]);
ylim([-0.1 0.8]);

set(gca,'fontsize', 14, 'fontname', 'times');

% Save figure
print(gcf, '../Figures/Figure_2.png', '-dpng', '-r300');

%% Figure 2 Replicate
% Define Model Parameters
p.alpha = 1; p.beta = 1;
p.a = 10; p.b = -10;
p.c = 10; p.d = 2;

% Find u* and v* such that θu = -2 & θv = -4
p.theta_u = -2;
p.theta_v = -4;
[p.u, p.v] = calcBias(p);

% Specify ranges for tau1, tau2, and omega
p.tau1_vals = linspace(0, 12, 121); % Range of tau1
p.tau2_vals = linspace(0, 1.5, 3); % Range of tau2
p.omega_vals = linspace(0, 2.5, 7); % Omega range

% Calculate hopf bifurcation lines in \tau_{1} and \tau_{2}
bifn = ddeBifn(p);

% Initialise figure 3
figure(3);
clf; hold on;

% Define a custom colour palette
% p.c_palette = customColourPalette(["f77e1b"; "2c70b3"; "c74440"]);
p.c_palette = customColourPalette(["000000"; "f77e1b"]);

% Define omega range
omega.min = min(cellfun(@min, bifn.omegas));
omega.max = max(cellfun(@max, bifn.omegas));

% Plot first bifurcation line
for i = 1:length(bifn.line_1) - 1
    % Define a segment
    x_seg = [bifn.line_1(1,i), bifn.line_1(1,i+1)];
    y_seg = [bifn.line_1(2,i), bifn.line_1(2,i+1)];
    
    % Map omega to custom colour pallete
    omega.val = bifn.line_1a(i+1);
    omega.colour_idx = round( ...
        ((omega.val - omega.min) / (omega.max - omega.min)) * (255)) + 1;
    omega.colour = p.c_palette(omega.colour_idx, :);
    
    % Plot the segment with the corresponding colour
    plot(x_seg, y_seg, 'color', omega.colour, 'linewidth', 2);
end
clear i x_seg y_seg

% Plot second bifurcation line
for i = 1:length(bifn.line_2)-1
    % Define a segment
    x_seg = [bifn.line_2(1,i), bifn.line_2(1,i+1)];
    y_seg = [bifn.line_2(2,i), bifn.line_2(2,i+1)];
    
    % Map omega to custom colour pallete
    omega.val = bifn.line_2a(i+1);
    omega.colour_idx = round( ...
        ((omega.val - omega.min) / (omega.max - omega.min)) * (255)) + 1;
    omega.colour = p.c_palette(omega.colour_idx, :);
    
    % Plot the segment with the corresponding colour
    plot(x_seg, y_seg, 'color', omega.colour, 'linewidth', 2);
end
clear i x_seg y_seg

% Plot label
plot(0.4, 1.44, '*', 'color', 'k', 'markersize', 12);
plot(0.4, 1.44, 'square', 'color', 'k', 'markersize', 20);

% Plot simulation points
plot(0.5, 1, '*', 'color', '#c74440', 'markersize', 5);
plot(3, 1, '*', 'color', '#2c70b3', 'markersize', 5);
plot(6, 1, '*', 'color', '#378c47', 'markersize', 5);
plot(10, 1, '*', 'color', '#f77e1b', 'markersize', 5);

% Format axes
xlabel("\tau_1");
ylabel("\tau_2", 'rotation', 0);

xlim([0 12]);
xticks([0 4 8 12]);
xticklabels({'0', '4', '8', '12'});

ylim([0 1.5]);
yticks([0 0.5 1.0 1.5]);
yticklabels({'0', '0.5', '1.0', '1.5'});

set(gca,'fontsize', 14, 'fontname', 'times');

% Add a colour bar for reference
colormap(p.c_palette);
clim([omega.min omega.max]);
colorbar;
clear omega

% Add annotations
annotation('textbox',[0.1696 0.4562 0.1219 0.0619],'string',"unstable", ...
    'rotation',71.5,'fontsize',14,'fontname','times','edgecolor','none');
annotation('textbox',[0.2428 0.4062 0.0969 0.0619],'string',"stable", ...
    'rotation',71.5,'fontsize',14,'fontname','times','edgecolor','none');
annotation('textbox',[0.4781 0.3905 0.1219 0.0619],'string',"unstable", ...
    'fontsize',14,'fontname','times','edgecolor','none');

% Save figure
print(gcf, '../Figures/Figure_3.png', '-dpng', '-r300');

%% Simulate DDE System around Bifurcation Lines
% Define Model Parameters
p.alpha = 1; p.beta = 1;
p.a = 10; p.b = -10;
p.c = 10; p.d = 2;
p.theta_u = -2; p.theta_v = -4;

% Define DDE parameters
p.tspan = [0 60];
p.history = [0.6 0.7];
p.options = ddeset('reltol', 1e-5);

p.tau1 = 0.5; p.tau2 = 1;
p.delays = [p.tau1 p.tau2];
[sol.dde_sim_1, nullclines] = ddeSim(p);

p.tau1 = 3; p.tau2 = 1;
p.delays = [p.tau1 p.tau2];
sol.dde_sim_2 = ddeSim(p);

p.tau1 = 6; p.tau2 = 1;
p.delays = [p.tau1 p.tau2];
sol.dde_sim_3 = ddeSim(p);

p.tau1 = 10; p.tau2 = 1;
p.delays = [p.tau1 p.tau2];
sol.dde_sim_4 = ddeSim(p);

% Initialise figure 4
figure(4); clf;
tiledlayout(2, 2, 'tilespacing', 'compact', 'padding', 'loose');

% First subplot
nexttile; hold on;

% Plot nullclines
plot(nullclines.u_null, nullclines.v_range, ...
    '-', 'color', '#c74440', 'linewidth', 1.5);
plot(nullclines.u_range, nullclines.v_null, ...
    '-', 'color', '#2c70b3', 'linewidth', 1.5);

% Plot solution
plot(sol.dde_sim_1.y(1,1), sol.dde_sim_1.y(2,1), ...
    '.', 'color', '#378c47', 'markersize', 8);
plot(sol.dde_sim_1.y(1,:), sol.dde_sim_1.y(2,:), ...
    'color', '#378c47', 'linewidth', 1);

% Plot label
plot(0.95, 0.94, '*', 'color', '#c74440', 'markersize', 12);
plot(0.95, 0.94, 'square', 'color', '#c74440', 'markersize', 20);

% Format axes
ylabel("$\mathit{v}$", 'interpreter', 'latex', 'rotation', 0, ...
    'position', [-0.2, 0.5, -1], 'verticalalignment', 'middle');

xlim([0 1]);
xticks(0:0.2:1);
xticklabels({'0','0.2','0.4','0.6','0.8','1.0'});

ylim([0 1]);
yticks(0:0.2:1);
yticklabels({'0','0.2','0.4','0.6','0.8','1.0'});

set(gca,'fontsize', 14, 'fontname', 'times');

% Second subplot
nexttile; hold on;

% Plot nullclines
plot(nullclines.u_null, nullclines.v_range, ...
    '-', 'color', '#c74440', 'linewidth', 1.5);
plot(nullclines.u_range, nullclines.v_null, ...
    '-', 'color', '#2c70b3', 'linewidth', 1.5);

% Plot solution
plot(sol.dde_sim_2.y(1,1), sol.dde_sim_2.y(2,1), ...
    '.', 'color', '#378c47', 'markersize', 8);
plot(sol.dde_sim_2.y(1,:), sol.dde_sim_2.y(2,:), ...
    'color', '#378c47', 'linewidth', 1);

% Plot label
plot(0.95, 0.94, '*', 'color', '#2c70b3', 'markersize', 12);
plot(0.95, 0.94, 'square', 'color', '#2c70b3', 'markersize', 20);

% Format axes
xlim([0 1]);
xticks(0:0.2:1);
xticklabels({'0','0.2','0.4','0.6','0.8','1.0'});

ylim([0 1]);
yticks(0:0.2:1);
yticklabels({'0','0.2','0.4','0.6','0.8','1.0'});

set(gca,'fontsize', 14, 'fontname', 'times');

% Third subplot
nexttile; hold on;

% Plot nullclines
plot(nullclines.u_null, nullclines.v_range, ...
    '-', 'color', '#c74440', 'linewidth', 1.5);
plot(nullclines.u_range, nullclines.v_null, ...
    '-', 'color', '#2c70b3', 'linewidth', 1.5);

% Plot solution
plot(sol.dde_sim_3.y(1,1), sol.dde_sim_3.y(2,1), ...
    '.', 'color', '#378c47', 'markersize', 8);
plot(sol.dde_sim_3.y(1,:), sol.dde_sim_3.y(2,:), ...
    'color', '#378c47', 'linewidth', 1);

% Plot label
plot(0.95, 0.94, '*', 'color', '#378c47', 'markersize', 12);
plot(0.95, 0.94, 'square', 'color', '#378c47', 'markersize', 20);

% Format axes
xlabel("$\mathit{u}$", 'interpreter', 'latex');
ylabel("$\mathit{v}$", 'interpreter', 'latex', 'rotation', 0, ...
    'position', [-0.2, 0.5, -1], 'verticalalignment', 'middle');

xlim([0 1]);
xticks(0:0.2:1);
xticklabels({'0','0.2','0.4','0.6','0.8','1.0'});

ylim([0 1]);
yticks(0:0.2:1);
yticklabels({'0','0.2','0.4','0.6','0.8','1.0'});

set(gca,'fontsize', 14, 'fontname', 'times');

% Fourth subplot
nexttile; hold on;

% Plot nullclines
plot(nullclines.u_null, nullclines.v_range, ...
    '-', 'color', '#c74440', 'linewidth', 1.5);
plot(nullclines.u_range, nullclines.v_null, ...
    '-', 'color', '#2c70b3', 'linewidth', 1.5);

% Plot solution
plot(sol.dde_sim_4.y(1,1), sol.dde_sim_4.y(2,1), ...
    '.', 'color', '#378c47', 'markersize', 8);
plot(sol.dde_sim_4.y(1,:), sol.dde_sim_4.y(2,:), ...
    'color', '#378c47', 'linewidth', 1);

% Plot label
plot(0.95, 0.94, '*', 'color', '#f77e1b', 'markersize', 12);
plot(0.95, 0.94, 'square', 'color', '#f77e1b', 'markersize', 20);

% Format axes
xlabel("$\mathit{u}$", 'interpreter', 'latex');

xlim([0 1]);
xticks(0:0.2:1);
xticklabels({'0','0.2','0.4','0.6','0.8','1.0'});

ylim([0 1]);
yticks(0:0.2:1);
yticklabels({'0','0.2','0.4','0.6','0.8','1.0'});

set(gca,'fontsize', 14, 'fontname', 'times');

% Save figure
print(gcf, '../Figures/Figure_4.png', '-dpng', '-r300');

%% Figure 3 Replicate
% Define Model Parameters
p.alpha = 1; p.beta = 600; % high \beta approximates heaviside
p.a = -1; p.b = -0.4;
p.c = -0.4; p.d = -1;
p.theta_u = 0.7; p.theta_v = 0.7;
p.tau1 = 1; p.tau2 = 1.4;

% Define DDE parameters
p.tspan = [0 15];
p.delays = [p.tau1 p.tau2];
p.options = ddeset('reltol', 1e-5);

% Solve DDE
% sol.Synchronous = dde23(@(t, y, Z) ddefun(t, y, Z, p), ...
%     p.delays,p.history1,p.tspan,p.options);
% sol.AntiSynchronous = dde23(@(t, y, Z) ddefun(t, y, Z, p), ...
%     p.delays,p.history2,p.tspan,p.options);

p.history = [0.37 0.8];
sol.Synchronous = ddeSim(p);

p.history = [0.2 0.8];
sol.AntiSynchronous = ddeSim(p);

% Initilaise figure 5
figure(5); clf;
tiledlayout(2,1,'tilespacing','compact','padding','loose');

% First subplot
nexttile; hold on;

% Plot synchronous solution
plot(sol.Synchronous.x, sol.Synchronous.y(1,:), ...
    'k', 'linewidth', 1.5); % u against time
dashline(sol.Synchronous.x,sol.Synchronous.y(2,:), ...
    1.5, 1, 1.5, 1, 'color', '#808080', 'linewidth', 1.5); % v against time

% Format axes
xlim([0, 15]);
xticks([0 5 10 15]);
xticklabels({'','','',''});

ylim([0, 1]);
yticks([0 0.5 1.0]);
yticklabels({'0', '0.5', '1.0'});

set(gca,'fontsize', 14, 'fontname', 'times');

% Second subplot
nexttile; hold on;

% Plot antisynchronous solution
plot(sol.AntiSynchronous.x, sol.AntiSynchronous.y(1,:), ...
    'k', 'linewidth', 1.5); % u against time
dashline(sol.AntiSynchronous.x,sol.AntiSynchronous.y(2,:), ...
    1.5, 1, 1.5, 1, 'color', '#808080', 'linewidth', 1.5); % v against time

% Format axes
xlabel("$\mathit{t}$", 'interpreter', 'latex');

xlim([0, 15]);
xticks([0 5 10 15]);
xticklabels({'0', '5', '10', '15'});

ylim([0, 1]);
yticks([0 0.5 1.0]);
yticklabels({'0', '0.5', '1.0'});

set(gca,'fontsize', 14, 'fontname', 'times');

% Add annotations
annotation('textbox',[0.14 0.81 0.0446 0.0536],'string',"$\mathit{u}$", ...
    'fontsize',14,'edgecolor','none','interpreter','latex');
annotation('textbox',[0.15 0.61 0.0446 0.0536],'String',"$\mathit{v}$", ...
    'fontsize',14,'edgecolor','none','interpreter','latex');
annotation('textbox',[0.15 0.37 0.0446 0.0536],'string',"$\mathit{u}$", ...
    'fontsize',14,'edgecolor','none','interpreter','latex');
annotation('textbox',[0.15 0.18 0.0446 0.0536],'string',"$\mathit{v}$", ...
    'fontsize',14,'edgecolor','none','interpreter','latex');

% Save figure
print(gcf, '../Figures/Figure_5.png', '-dpng', '-r300');

%% Figure 6 Replicate
% Define model parameters
p.alpha = 1; p.beta = 60;
p.a = -1; p.b = -0.4;
p.c = -1; p.d = 0;
p.theta_u = 0.7; p.theta_v = 0.5;
p.tau_1 = 0.2; p.tau_2 = p.tau_1;

bifn2 = ddeBiftoolMain(p);

% Initialise figure 6
figure(6);
clf; hold on;

% Plot bifurcations
    % First saddle-node bifurcation
plot(bifn2.sn1.x, bifn2.sn1.y, 'k', 'linewidth', 1);
    % Second saddle-node bifurcation
plot(bifn2.sn2.x, bifn2.sn2.y, 'k', 'linewidth', 1);
    % Hopf bifurcation
plot(bifn2.hopf.x, bifn2.hopf.y, 'k--', 'linewidth', 1);
    % First saddle-node bifurcation of periodic orbits
plot(bifn2.snpo1.x, bifn2.snpo1.y, 'k-o', 'markersize', 4, 'linewidth', 1);
    % Second saddle-node bifurcation of periodic orbits
plot(bifn2.snpo2.x, bifn2.snpo2.y, 'k-o', 'markersize', 4, 'linewidth', 1);

% Format axes
xlabel("$\theta_{\mathit{u}}$", 'interpreter', 'latex');
ylabel('\tau', 'rotation', 0, 'position', [0.34, 0.25, -1], ...
    'verticalalignment', 'middle');

xlim([0.4, 1]);
xticks([0.4 0.5 0.6 0.7 0.8 0.9 1]);
xticklabels({'0.4','0.5','0.6','0.7', '0.8', '0.9', '1.0'});

ylim([0, 0.5]);
yticks([0 0.1 0.2 0.3 0.4 0.5]);
yticklabels({'0', '0.1', '0.2', '0.3', '0.5', '0.5'});

set(gca,'fontsize', 14, 'fontname', 'times');

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
figure(7); clf;
tiledlayout(3,1,'tilespacing','compact','padding','loose');

% Subplot 1 - \tau = 0.5
nexttile; hold on;

plot(stst.a.stable1.x, stst.a.stable1.y, 'k-');
plot(stst.a.unstable.x, stst.a.unstable.y, 'k--');
plot(stst.a.stable2.x, stst.a.stable2.y, 'k-');

if isfield(po.a, 'stable1')
    plot(po.a.stable1.x, po.a.stable1.y, 'k-o', 'markersize', 4);
end
if isfield(po.a, 'unstable')
    plot(po.a.unstable.x, po.a.unstable.y, 'k-x', 'markersize', 4);
end
if isfield(po.a, 'stable2')
    plot(po.a.stable2.x, po.a.stable2.y, 'k-o', 'markersize', 4);
end

% Format axes
ylabel("$\mathit{u}$", 'rotation', 0, 'interpreter', 'latex', ...
    'position', [0.35, 0.5, -1], 'verticalalignment', 'middle');

xlim([0.4 1]);
xticks([0.4 0.5 0.6 0.7 0.8 0.9 1]);
xticklabels({'0.4','0.5','0.6','0.7', '0.8', '0.9', '1.0'});

ylim([0 1]);
yticks([0 0.5 1]);
yticklabels({'0', '0.5', '1.0'});

set(gca,'fontsize', 14, 'fontname', 'times');

% Subplot 2 - \tau = 0.2
nexttile; hold on;

plot(stst.b.stable1.x, stst.b.stable1.y, 'k-');
plot(stst.b.unstable.x, stst.b.unstable.y, 'k--');
plot(stst.b.stable2.x, stst.b.stable2.y, 'k-');

if isfield(po.b, 'stable1')
    plot(po.b.stable1.x, po.b.stable1.y, 'k-o', 'markersize', 4);
end
if isfield(po.b, 'unstable')
    plot(po.b.unstable.x, po.b.unstable.y, 'k-x', 'markersize', 4);
end
if isfield(po.b, 'stable2')
    plot(po.b.stable2.x, po.b.stable2.y, 'k-o', 'markersize', 4);
end

% Format axes
ylabel("$\mathit{u}$", 'rotation', 0, 'interpreter', 'latex', ...
    'position', [0.465, 0.45, -1], 'verticalalignment', 'middle');

xlim([0.5 0.9]);
xticks([0.5 0.6 0.7 0.8 0.9]);
xticklabels({'0.5','0.6','0.7','0.8', '0.9'});

ylim([0 0.9]);
yticks([0.2 0.4 0.6 0.8]);
yticklabels({'0.2', '0.4', '0.6', '0.8'});

set(gca,'fontsize', 14, 'fontname', 'times');

% Subplot 3 - \tau = 0.09
nexttile; hold on;

plot(stst.c.stable1.x, stst.c.stable1.y, 'k-');
plot(stst.c.unstable.x, stst.c.unstable.y, 'k--');
plot(stst.c.stable2.x, stst.c.stable2.y, 'k-');

if isfield(po.c, 'stable1')
    plot(po.c.stable1.x, po.c.stable1.y, 'k-o', 'markersize', 4);
end
if isfield(po.c, 'unstable')
    plot(po.c.unstable.x, po.c.unstable.y, 'k-x', 'markersize', 4);
end
if isfield(po.c, 'stable2')
    plot(po.c.stable2.x, po.c.stable2.y, 'k-o', 'markersize', 4);
end

% Format axes
xlabel("$\theta_{\mathit{u}}$", 'interpreter', 'latex');
ylabel("$\mathit{u}$", 'rotation', 0, 'interpreter', 'latex', ...
    'position', [0.465, 0.5, -1], 'verticalalignment', 'middle');

xlim([0.5 0.9]);
xticks([0.5 0.6 0.7 0.8 0.9]);
xticklabels({'0.5','0.6','0.7', '0.8', '0.9'});

ylim([0.3 0.7]);
yticks([0.3 0.4 0.5 0.6 0.7]);
yticklabels({'0.3', '0.4', '0.5', '0.6', '0.7'});

set(gca,'fontsize', 14, 'fontname', 'times');

% Save figure
print(gcf, '../Figures/Figure_7.png', '-dpng', '-r300');

%% Figure 9a Replicate
% Select Parameters
p.alpha = 1; p.beta = 60;
p.a = -6; p.b = 2.5;
p.c = 2.5; p.d = -6;
p.theta_u = 0.2; p.theta_v = 0.2;
p.tau1 = 0.1; p.tau2 = 0.1;

% Define DDE parameters
p.tspan = [0 40];
p.delays = [p.tau1 p.tau2];
p.history = [0.074 0.077];
p.options = ddeset('reltol', 1e-5);

% Solve DDE
sol.Chaotic = ddeSim(p);

% Initialise figure 8
figure(8);
clf; hold on;

% Plot chotic solution
plot(sol.Chaotic.y(1,:), sol.Chaotic.y(2,:), ...
    'k', 'linewidth', 1.5); % u against v

% Format axes
xlabel("$\mathit{u}$", 'interpreter', 'latex');
ylabel("$\mathit{v}$", 'interpreter', 'latex','rotation', 0);

xlim([0.062, 0.178]);
xticks([0.08 0.1 0.12 0.14 0.16]);
yticks([0.06 0.08 0.1 0.12 0.14 0.16 0.18]);

ylim([0.06, 0.18]);
xticklabels({'0.08', '', '0.12', '', '0.16', ''});
yticklabels({'0.06', '', '0.10', '', '0.14', '', '0.18'});

set(gca,'fontsize', 14, 'fontname', 'times');

% Save figure
print(gcf, '../Figures/Figure_9a.png', '-dpng', '-r300');

%% Figure 9b Replicate
% Select Parameters
p.alpha = 1; p.beta = 40;
p.a = -6; p.b = 2.5;
p.c = 2.5; p.d = -6;
p.theta_u = 0.2; p.theta_v = 0.2;
p.tau1 = 0.1; p.tau2 = 0.1;

% Define DDE parameters
p.tspan = [0 40];
p.delays = [p.tau1 p.tau2];
p.history = [0.074 0.077];
p.options = ddeset('reltol', 1e-5);

% Solve DDE
sol.QuasiPeriodic = ddeSim(p);

% Initialise figure 9
figure(9);
clf; hold on;

% Plot chotic solution
plot(sol.QuasiPeriodic.y(1,:), sol.QuasiPeriodic.y(2,:), ...
    'k', 'linewidth', 1.5); % u against v

% Format axes
xlabel("$\mathit{u}$", 'interpreter', 'latex');
ylabel("$\mathit{v}$", 'interpreter', 'latex','rotation', 0);

xlim([0.07, 0.105]);
xticks([0.07 0.075 0.08 0.085 0.09 0.095 0.1 0.105]);
xticklabels({'', '', '0.08', '', '0.09', '', '0.10', ''});

ylim([0.07, 0.105]);
yticks([0.07 0.075 0.08 0.085 0.09 0.095 0.1 0.105]);
yticklabels({'', '0.075', '', '0.085', '', '0.095', '', '0.105'});

set(gca,'fontsize', 14, 'fontname', 'times');

% Save figure
print(gcf, '../Figures/Figure_9b.png', '-dpng', '-r300');

%% Calculate Maximum Lyapunov Exponents over parameter space
% This section may take up to 3 hours to complete

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
p.options = ddeset('reltol', 1e-5);

% Calculate maximum Lyapunov Exponents
% lambda_max = calcLyapunovExponent(p);

% Alternatively, load pre-calculated lambda_max
load('lambda_max.mat');
%% Figure 8 Replictae
% Extract and filter lambda_max
le_img = lambda_max;
le_img(le_img < 0) = 0;

% Initialise figure 10
figure(10);
clf; hold on;

% Plot lambda max
imagesc(p.a_grid(1, :), p.b_grid(:, 1), le_img);
colormap(flipud(gray));
colorbar;

% Calculate plotting boundaries
astep = (p.arange(2) - p.arange(1)) / 2;
bstep = (p.brange(2) - p.brange(1)) / 2;

% Format axes
xlabel("$\mathit{a}$", 'interpreter', 'latex');
ylabel("$\mathit{b}$", 'interpreter', 'latex','rotation', 0);

xlim([-10 - astep, 0 + astep]);
xticks([-10 -8 -6 -4 -2 0]);
xticklabels({'-10', '-8', '-6', '-4', '-2', '0'});

ylim([0 - bstep, 5 + bstep]);
yticks([0 1 2 3 4 5]);
yticklabels({'0', '1', '2', '3', '4', '5'});

set(gca,'fontsize', 14, 'fontname', 'times');

% Save figure
print(gcf, '../Figures/Figure_8.png', '-dpng', '-r300');

clear astep bstep