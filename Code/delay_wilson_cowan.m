%% Modelling Delays in an Activity-Based Wilson-Cowan Neural Network
clearvars; clc;
addpath('helper_functions/');

%% Figure 1 Replicate - With Bogdanov-Takens Bifurcations
% Define Model Parameters
p.alpha = 1; p.beta = 1;
p.a = 10; p.b = -10;
p.c = 10; p.d = 2;

% Calculate hopf, saddle-node and bogdanov-takens bifurcations
[hopf, sn, bt] = odeBifn(p);

% Initialise figure 1
figure(1);
set(gcf, 'position', [400, 100, 600, 450]);
clf; hold on;

% Plot hopf bifurcation
plot(hopf.theta.uP, hopf.theta.vP, 'k--', 'linewidth',2); % plus values
plot(hopf.theta.uM, hopf.theta.vM,'k--', 'linewidth',2); % minus values

% Plot saddle node bifurcation
plot(sn.theta.uP, sn.theta.vP, 'k', 'linewidth',3); % plus values
plot(sn.theta.uM, sn.theta.vM, 'k', 'linewidth',3); % minus values

% Plot bogdanov-takens bifurcations
plot(bt.theta.uP, bt.theta.vP, 'o', 'color', '#c74440', ...
    'markersize', 16, 'linewidth', 2); % plus values
plot(bt.theta.uM, bt.theta.vM, 'o', 'color', '#c74440', ...
    'markersize', 16, 'linewidth', 2); % minus values

% Plot simulation points
% plot(-4.1, -7, '*', 'color', '#c74440', 'markersize', 5);
% plot(-3.7, -7, '*', 'color', '#2c70b3', 'markersize', 5);
% plot(-3.3, -7, '*', 'color', '#378c47', 'markersize', 5);
% plot(-2.9, -7, '*', 'color', '#f77e1b', 'markersize', 5);

% Plot figure 3 point
plot(-2, -4, '*', 'color', '#f77e1b', 'markersize', 16, 'linewidth', 2);
text(-3.1, -3.3, "fig. 3a", 'fontsize', 20, ...
    'fontname', 'times', 'FontAngle', 'italic');

% Define curves bounding oscillatory region
boundary.u1 = hopf.theta.uP(1:761);
boundary.v1 = hopf.theta.vP(1:761);
boundary.u2 = hopf.theta.uP(11579:end);
boundary.v2 = hopf.theta.vP(11579:end);
boundary.u3 = hopf.theta.uM(1:6946);
boundary.v3 = hopf.theta.vM(1:6946);
boundary.u4 = hopf.theta.uM(17764:end);
boundary.v4 = hopf.theta.vM(17764:end);
boundary.u5 = sn.theta.uP(33:3186);
boundary.v5 = sn.theta.vP(33:3186);
boundary.u6 = sn.theta.uM(74274:77427);
boundary.v6 = sn.theta.vM(74274:77427);

% Concatenate boundaries into one curve
boundary.u = [fliplr(boundary.u1), boundary.u3, ...
     fliplr(boundary.u5), boundary.u4, ...
     fliplr(boundary.u2), boundary.u6];
boundary.v = [fliplr(boundary.v1), boundary.v3, ...
     fliplr(boundary.v5), boundary.v4, ...
     fliplr(boundary.v2), boundary.v6];

% Plot oscillatory region
fill(boundary.u, boundary.v, [0.8, 0.8, 0.8], ...
    'facealpha', 0.3, 'edgecolor', 'none'); % Light grey shading

% Format axes
xlabel("$\theta_{\mathit{u}}$", 'interpreter', 'latex');
ylabel("$\theta_{\mathit{v}}$", 'rotation', 0, 'interpreter', 'latex');

xlim([-7 7]);
xticks(-6:2:6);
xticklabels({'','-4','','0','','4',''});

ylim([-12 0]);
yticks(-12:2:0);
yticklabels({'-12','','-8','','-4','','0'});

set(gca,'fontsize', 24, 'fontname', 'times');

% Graph Text
text(0, -9.2, "HB", 'fontsize', 24, 'fontname', 'times');
text(5, -4, "SN", 'fontsize', 24, 'fontname', 'times');
text(-0.5, -6, "OR", 'fontsize', 24, 'fontname', 'times');

% Save figure
print(gcf, '../Figures/Figure_1.png', '-dpng', '-r300');

%% Simulate ODE System Around Bifurcation Lines
% Define Model Parameters
p.alpha = 1; p.beta = 1;
p.a = 10; p.b = -10;
p.c = 10; p.d = 2;

% Define simulation parameters
p.tspan = [0 60]; p.history = [0.7 0.1];

% Initialise figure 2
figure(2);
set(gcf, 'position', [400, 100, 600, 450]);
clf;
tiledlayout(2,2,'tilespacing','compact','padding','loose');

% Subplot 1
nexttile; hold on;

% Calculate nullclines and simulate
p.theta_u = -4.1; p.theta_v = -7;
[~, ~, sol.ode_sim_1, nullclines] = odeSim(p);

% Plot nullclines
plot(nullclines.u_null, nullclines.v_range, ...
    '-.', 'color', '#c74440', 'linewidth', 3);
plot(nullclines.u_range, nullclines.v_null, ...
    '-.', 'color', '#2c70b3', 'linewidth', 3);

% Plot simulation
plot(sol.ode_sim_1(1,1), sol.ode_sim_1(1,2), ...
    '.', 'color', '#378c47', 'markersize', 12);
plot(sol.ode_sim_1(:,1), sol.ode_sim_1(:,2), ...
    '-', 'color', '#378c47', 'linewidth', 2);
plot(sol.ode_sim_1(end,1), sol.ode_sim_1(end,2), ...
    '.', 'color', 'k', 'markersize', 16);

% Plot label
annotation('textbox',[0.132,0.895,0.06,0.06],'string',"(a)", ...
    'fontsize',24,'fontname','times','edgecolor','none');
plot(0.75, 0.73, '*', 'color', '#c74440', ...
    'markersize', 16, 'linewidth', 2);
plot(0.75, 0.73, 'square', 'color', '#c74440', ...
    'markersize', 24, 'linewidth', 2);

% Format axes
ylabel("$\mathit{v}$", 'interpreter', 'latex', 'rotation', 0, ...
    'position', [-0.32, 0.35, -1], 'verticalalignment', 'middle');

xlim([-0.1 0.81]);
ylim([-0.1 0.81]);

set(gca,'fontsize', 24, 'fontname', 'times');

% Subplot 2
nexttile; hold on;

% Calculate nullclines and simulate
p.theta_u = -3.7; p.theta_v = -7;
[~, ~, sol.ode_sim_2, nullclines] = odeSim(p);

% Plot nullclines
plot(nullclines.u_null, nullclines.v_range, ...
    '-.', 'color', '#c74440', 'linewidth', 3);
plot(nullclines.u_range, nullclines.v_null, ...
    '-.', 'color', '#2c70b3', 'linewidth', 3);

% Plot simulation
plot(sol.ode_sim_2(1,1), sol.ode_sim_2(1,2), ...
    '.', 'color', '#378c47', 'markersize', 12);
plot(sol.ode_sim_2(:,1), sol.ode_sim_2(:,2), ...
    '-', 'color', '#378c47', 'linewidth', 2);
plot(sol.ode_sim_2(end,1), sol.ode_sim_2(end,2), ...
    '.', 'color', 'k', 'markersize', 16);

% Plot label
annotation('textbox',[0.553,0.895,0.06,0.06],'string',"(b)", ...
    'fontsize',24,'fontname','times','edgecolor','none');
plot(0.75, 0.73, '*', 'color', '#2c70b3', ...
    'markersize', 16, 'linewidth', 2);
plot(0.75, 0.73, 'square', 'color', '#2c70b3', ...
    'markersize', 24, 'linewidth', 2);

% Format axes
xlim([-0.1 0.81]);
ylim([-0.1 0.81]);

set(gca,'fontsize', 24, 'fontname', 'times');

% Subplot 3
nexttile; hold on;

% Calculate nullclines and simulate
p.theta_u = -3.3; p.theta_v = -7;
[~, ~, sol.ode_sim_3, nullclines] = odeSim(p);

% Plot nullclines
plot(nullclines.u_null, nullclines.v_range, ...
    '-.', 'color', '#c74440', 'linewidth', 3);
plot(nullclines.u_range, nullclines.v_null, ...
    '-.', 'color', '#2c70b3', 'linewidth', 3);

% Plot simulation
plot(sol.ode_sim_3(1,1), sol.ode_sim_3(1,2), ...
    '.', 'color', '#378c47', 'markersize', 12);
plot(sol.ode_sim_3(:,1), sol.ode_sim_3(:,2), ...
    '-', 'color', '#378c47', 'linewidth', 2);
plot(sol.ode_sim_3(end,1), sol.ode_sim_3(end,2), ...
    '.', 'color', 'k', 'markersize', 16);

% Plot label
annotation('textbox',[0.132,0.450,0.06,0.06],'string',"(c)", ...
    'fontsize',24,'fontname','times','edgecolor','none');
plot(0.75, 0.73, '*', 'color', '#378c47', ...
    'markersize', 16, 'linewidth', 2);
plot(0.75, 0.73, 'square', 'color', '#378c47', ...
    'markersize', 24, 'linewidth', 2);

% Format axes
xlabel("$\mathit{u}$", 'interpreter', 'latex');
ylabel("$\mathit{v}$", 'interpreter', 'latex', 'rotation', 0, ...
    'position', [-0.32, 0.35, -1], 'verticalalignment', 'middle');

xlim([-0.1 0.81]);
ylim([-0.1 0.81]);

set(gca,'fontsize', 24, 'fontname', 'times');

% Subplot 4
nexttile; hold on;

% Calculate nullclines and simulate
p.theta_u = -2.9; p.theta_v = -7;
[~, ~, sol.ode_sim_4, nullclines] = odeSim(p);

% Plot nullclines
plot(nullclines.u_null, nullclines.v_range, ...
    '-.', 'color', '#c74440', 'linewidth', 3);
plot(nullclines.u_range, nullclines.v_null, ...
    '-.', 'color', '#2c70b3', 'linewidth', 3);

% Plot simulation
plot(sol.ode_sim_4(1,1), sol.ode_sim_4(1,2), ...
    '.', 'color', '#378c47', 'markersize', 12);
plot(sol.ode_sim_4(:,1), sol.ode_sim_4(:,2), ...
    '-', 'color', '#378c47', 'linewidth', 2);

% Plot label
annotation('textbox',[0.553,0.450,0.06,0.06],'string',"(d)", ...
    'fontsize',24,'fontname','times','edgecolor','none');
plot(0.75, 0.73, '*', 'color', '#f77e1b', ...
    'markersize', 16, 'linewidth', 2);
plot(0.75, 0.73, 'square', 'color', '#f77e1b', ...
    'markersize', 24, 'linewidth', 2);

% Format axes
xlabel("$\mathit{u}$", 'interpreter', 'latex');

xlim([-0.1 0.81]);
ylim([-0.1 0.81]);

set(gca,'fontsize', 24, 'fontname', 'times');

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
set(gcf, 'position', [400, 100, 600, 300]);
clf; hold on;

% Define a custom colour palette
p.c_palette = customColourPalette(["000000"; "f77e1b"]);

% Define omega range
omega.min = 0;
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
    plot(x_seg, y_seg, 'color', omega.colour, 'linewidth', 4);
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
    plot(x_seg, y_seg, 'color', omega.colour, 'linewidth', 4);
end
clear i x_seg y_seg

% Plot simulation points
plot(0.5, 1, '*', 'color', '#c74440', 'markersize', 16, 'linewidth', 2);
text(0.2, 1.17, "fig. 3b", 'fontsize', 20, ...
    'fontname', 'times', 'FontAngle', 'italic');

plot(3, 1, '*', 'color', '#2c70b3', 'markersize', 16, 'linewidth', 2);
text(2.28, 1.17, "fig. 3c", 'fontsize', 20, ...
    'fontname', 'times', 'FontAngle', 'italic');

plot(6, 1, '*', 'color', '#378c47', 'markersize', 16, 'linewidth', 2);
text(5, 1.17, "fig. 3d", 'fontsize', 20, ...
    'fontname', 'times', 'FontAngle', 'italic');


% Add a colour bar for reference
colormap(p.c_palette);
clim([0 2]);
colorbar;
set(colorbar, ...
    'Ticks', 0:0.5:2, ...
    'TickLabels', {'0.0', '0.5', '1.0', '1.5', '2.0'});
clear omega

% Format axes
xlabel("\tau_1");
ylabel("\tau_2", 'rotation', 0);

xlim([0 12]);
xticks([0 4 8 12]);
xticklabels({'0', '4', '8', '12'});

ylim([0 1.5]);
yticks([0 0.5 1.0 1.5]);
yticklabels({'0', '0.5', '1.0', '1.5'});

set(gca,'fontsize', 24, 'fontname', 'times');

% Add annotations
annotation('textbox',[0.1346 0.4262 0.1219 0.0619],'string',"unstable", ...
    'rotation',63.5,'fontsize',24,'fontname','times','edgecolor','none');
annotation('textbox',[0.2078 0.4062 0.0969 0.0619],'string',"stable", ...
    'rotation',63.5,'fontsize',24,'fontname','times','edgecolor','none');
annotation('textbox',[0.4081 0.4505 0.1219 0.0619],'string',"unstable", ...
    'fontsize',24,'fontname','times','edgecolor','none');

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

[~, ~, sol.ode_sim_5, nullclines] = odeSim(p);

p.tau_1 = 0.5; p.tau_2 = 1;
p.delays = [p.tau_1 p.tau_2];
[sol.dde_sim_1] = ddeSim(p);

p.tau_1 = 3; p.tau_2 = 1;
p.delays = [p.tau_1 p.tau_2];
sol.dde_sim_2 = ddeSim(p);

p.tau_1 = 6; p.tau_2 = 1;
p.delays = [p.tau_1 p.tau_2];
sol.dde_sim_3 = ddeSim(p);

% Initialise figure 4
figure(4);
set(gcf, 'position', [400, 100, 600, 700]);
clf;
tiledlayout(4, 2, 'tilespacing', 'compact', 'padding', 'loose');

% First subplot - (u,v) simulation 1
nexttile; hold on;

% Plot nullclines
plot(nullclines.u_null, nullclines.v_range, ...
    '-.', 'color', '#c74440', 'linewidth', 3);
plot(nullclines.u_range, nullclines.v_null, ...
    '-.', 'color', '#2c70b3', 'linewidth', 3);

% Plot solution
plot(sol.ode_sim_5(1,1), sol.ode_sim_5(1,2), ...
    '.', 'color', '#378c47', 'markersize', 16);
plot(sol.ode_sim_5(:,1), sol.ode_sim_5(:,2), ...
    '-', 'color', '#378c47', 'linewidth', 2);
plot(sol.ode_sim_5(end,1), sol.ode_sim_5(end,2), ...
    '.', 'color', 'k', 'markersize', 20);

% Plot label
annotation('textbox',[0.005,0.896,0.06,0.06],'string',"(a)", ...
    'fontsize',24,'fontname','times','edgecolor','none');
plot(0.95, 0.94, '*', 'color', '#f77e1b', ...
    'markersize', 20, 'linewidth', 2);
plot(0.95, 0.94, 'square', 'color', '#f77e1b', ...
    'markersize', 24, 'linewidth', 2);

% Format axes
xlabel("$\mathit{u}$", 'interpreter', 'latex');
ylabel("$\mathit{v}$", 'interpreter', 'latex', 'rotation', 0, ...
    'position', [-0.275, 0.5, -1], 'verticalalignment', 'middle');

xlim([0 1]);
xticks(0:0.5:1);
xticklabels({'0','0.5','1.0'});
xtickangle(0);

ylim([0 1]);
yticks(0:0.5:1);
yticklabels({'0','0.5','1.0'});

set(gca,'fontsize', 24, 'fontname', 'times');

% Second subplot - (u,v) simulation 2
nexttile; hold on;

% Plot nullclines
plot(nullclines.u_null, nullclines.v_range, ...
    '-.', 'color', '#c74440', 'linewidth', 3);
plot(nullclines.u_range, nullclines.v_null, ...
    '-.', 'color', '#2c70b3', 'linewidth', 3);

% Plot solution
plot(sol.dde_sim_1.y(1,1), sol.dde_sim_1.y(2,1), ...
    '.', 'color', '#378c47', 'markersize', 16);
plot(sol.dde_sim_1.y(1,:), sol.dde_sim_1.y(2,:), ...
    '-', 'color', '#378c47', 'linewidth', 2);

% Plot label
annotation('textbox',[0.93,0.896,0.06,0.06],'string',"(b)", ...
    'fontsize',24,'fontname','times','edgecolor','none');
plot(0.95, 0.94, '*', 'color', '#c74440', ...
    'markersize', 20, 'linewidth', 2);
plot(0.95, 0.94, 'square', 'color', '#c74440', ...
    'markersize', 24, 'linewidth', 2);

% Format axes
xlabel("$\mathit{u}$", 'interpreter', 'latex');
xlim([0 1]);
xticks(0:0.5:1);
xticklabels({'0','0.5','1.0'});
xtickangle(0);

ylim([0 1]);
yticks(0:0.5:1);
yticklabels({'0','0.5','1.0'});

set(gca,'fontsize', 24, 'fontname', 'times');

% Third subplot - u & v against time for simulation 1
nexttile; hold on;

sol.ode_sim_5_x = linspace(p.tspan(1), p.tspan(2), length(sol.ode_sim_5));

plot(sol.ode_sim_5_x, sol.ode_sim_5(:,1), ...
    '-', 'color', '#f77e1b', 'linewidth', 2);
plot(sol.ode_sim_5_x, sol.ode_sim_5(:,2), ...
    '-', 'color', '#6142a6', 'linewidth', 2);

% Format axes
xlabel("$\mathit{time, t}$", 'interpreter', 'latex');
ylabel("$\mathit{u,v}$", 'interpreter', 'latex', 'rotation', 90, ...
    'position', [-16, 0.5, -1], 'verticalalignment', 'middle')
xlim([0 60]);
xticks(0:15:60);
xticklabels({'0','','30','','60'});
xtickangle(0);

ylim([0 1]);
yticks(0:0.5:1);
yticklabels({'0','0.5','1.0'});

set(gca,'fontsize', 24, 'fontname', 'times');

% Fourth subplot - u & v againt time for simulation 2
nexttile; hold on;

plot(sol.dde_sim_1.x, sol.dde_sim_1.y(1,:), ...
    '-', 'color', '#f77e1b', 'linewidth', 2);
plot(sol.dde_sim_1.x, sol.dde_sim_1.y(2,:), ...
    '-', 'color', '#6142a6', 'linewidth', 2);

% Format axes
xlabel("$\mathit{time, t}$", 'interpreter', 'latex');
xlim([0 60]);
xticks(0:15:60);
xticklabels({'0','','30','','60'});
xtickangle(0);

ylim([0 1]);
yticks(0:0.5:1);
yticklabels({'0','0.5','1.0'});

set(gca,'fontsize', 24, 'fontname', 'times');

% Fifth subplot - (u,v) simulation 3
nexttile; hold on;

% Plot nullclines
plot(nullclines.u_null, nullclines.v_range, ...
    '-.', 'color', '#c74440', 'linewidth', 3);
plot(nullclines.u_range, nullclines.v_null, ...
    '-.', 'color', '#2c70b3', 'linewidth', 3);

% Plot solution
plot(sol.dde_sim_2.y(1,1), sol.dde_sim_2.y(2,1), ...
    '.', 'color', '#378c47', 'markersize', 16);
plot(sol.dde_sim_2.y(1,:), sol.dde_sim_2.y(2,:), ...
    '-', 'color', '#378c47', 'linewidth', 2);
plot(sol.dde_sim_2.y(1,end), sol.dde_sim_2.y(2,end), ...
    '.', 'color', '#000000', 'markersize', 20);

% Plot label
annotation('textbox',[0.005,0.426,0.06,0.06],'string',"(c)", ...
    'fontsize',24,'fontname','times','edgecolor','none');
plot(0.95, 0.94, '*', 'color', '#2c70b3', ...
    'markersize', 20, 'linewidth', 2);
plot(0.95, 0.94, 'square', 'color', '#2c70b3', ...
    'markersize', 24, 'linewidth', 2);

% Format axes
xlabel("$\mathit{u}$", 'interpreter', 'latex');
ylabel("$\mathit{v}$", 'interpreter', 'latex', 'rotation', 0, ...
    'position', [-0.275, 0.5, -1], 'verticalalignment', 'middle');

xlim([0 1]);
xticks(0:0.5:1);
xticklabels({'0','0.5','1.0'});
xtickangle(0);

ylim([0 1]);
yticks(0:0.5:1);
yticklabels({'0','0.5','1.0'});

set(gca,'fontsize', 24, 'fontname', 'times');

% Sixth subplot - (u,v) simulation 4
nexttile; hold on;

% Plot nullclines
plot(nullclines.u_null, nullclines.v_range, ...
    '-.', 'color', '#c74440', 'linewidth', 3);
plot(nullclines.u_range, nullclines.v_null, ...
    '-.', 'color', '#2c70b3', 'linewidth', 3);

% Plot solution
plot(sol.dde_sim_3.y(1,1), sol.dde_sim_3.y(2,1), ...
    '.', 'color', '#378c47', 'markersize', 16);
plot(sol.dde_sim_3.y(1,:), sol.dde_sim_3.y(2,:), ...
    '-', 'color', '#378c47', 'linewidth', 2);

% Plot label
annotation('textbox',[0.93,0.426,0.06,0.06],'string',"(d)", ...
    'fontsize',24,'fontname','times','edgecolor','none');
plot(0.95, 0.94, '*', 'color', '#378c47', ...
    'markersize', 20, 'linewidth', 2);
plot(0.95, 0.94, 'square', 'color', '#378c47', ...
    'markersize', 24, 'linewidth', 2);

% Format axes
xlabel("$\mathit{u}$", 'interpreter', 'latex');

xlim([0 1]);
xticks(0:0.5:1);
xticklabels({'0','0.5','1.0'});
xtickangle(0);

ylim([0 1]);
yticks(0:0.5:1);
yticklabels({'0','0.5','1.0'});

set(gca,'fontsize', 24, 'fontname', 'times');

% Seventh subplot - u & v against time for simulation 3
nexttile; hold on;

plot(sol.dde_sim_2.x, sol.dde_sim_2.y(1,:), ...
    '-', 'color', '#f77e1b', 'linewidth', 2);
plot(sol.dde_sim_2.x, sol.dde_sim_2.y(2,:), ...
    '-', 'color', '#6142a6', 'linewidth', 2);

% Format axes
xlabel("$\mathit{time, t}$", 'interpreter', 'latex');
ylabel("$\mathit{u,v}$", 'interpreter', 'latex', 'rotation', 90, ...
    'position', [-16, 0.5, -1], 'verticalalignment', 'middle')
xlim([0 60]);
xticks(0:15:60);
xticklabels({'0','','30','','60'});
xtickangle(0);

ylim([0 1]);
yticks(0:0.5:1);
yticklabels({'0','0.5','1.0'});

set(gca,'fontsize', 24, 'fontname', 'times');

% Eighth subplot - u & v againt time for simulation 4
nexttile; hold on;

plot(sol.dde_sim_3.x, sol.dde_sim_3.y(1,:), ...
    '-', 'color', '#f77e1b', 'linewidth', 2);
plot(sol.dde_sim_3.x, sol.dde_sim_3.y(2,:), ...
    '-', 'color', '#6142a6', 'linewidth', 2);

% Format axes
xlabel("$\mathit{time, t}$", 'interpreter', 'latex');
xlim([0 60]);
xticks(0:15:60);
xticklabels({'0','','30','','60'});
xtickangle(0);

ylim([0 1]);
yticks(0:0.5:1);
yticklabels({'0','0.5','1.0'});

set(gca,'fontsize', 24, 'fontname', 'times');

annotation('rectangle', [0.01, 0.48, 0.49, 0.47], ...
    'Color', 'k', 'LineWidth', 2);
annotation('rectangle', [0.5, 0.48, 0.49, 0.47], ...
    'Color', 'k', 'LineWidth', 2);
annotation('rectangle', [0.01, 0.01, 0.49, 0.47], ...
    'Color', 'k', 'LineWidth', 2);
annotation('rectangle', [0.5, 0.01, 0.49, 0.47], ...
    'Color', 'k', 'LineWidth', 2);

% Save figure
print(gcf, '../Figures/Figure_4.png', '-dpng', '-r300');

%% Figure 6 Replicate
% Define model parameters
p.alpha = 1; p.beta = 60;
p.a = -1; p.b = -0.4;
p.c = -1; p.d = 0;
p.theta_u = 0.7; p.theta_v = 0.5;
p.tau_1 = 0.2; p.tau_2 = p.tau_1;

bifn2 = ddeBiftoolMain(p);

% Initialise figure 5
figure(5);
set(gcf, 'position', [400, 100, 600, 350]);
clf; hold on;

% Plot bifurcations
    % First saddle-node bifurcation
plot(bifn2.sn1.x, bifn2.sn1.y, 'color', '#c74440', 'linewidth', 3);
    % Second saddle-node bifurcation
plot(bifn2.sn2.x, bifn2.sn2.y, 'color', '#c74440', 'linewidth', 3);
    % Hopf bifurcation
plot(bifn2.hopf.x, bifn2.hopf.y, '--', ...
    'color', '#2c70b3', 'linewidth', 3);
    % First saddle-node bifurcation of periodic orbits
plot(bifn2.snpo1.x, bifn2.snpo1.y, '-', ...
    'color', '#378c47', 'linewidth', 3);
plot(bifn2.snpo1.x, bifn2.snpo1.y, 'o', ...
    'color', '#378c47', 'markersize', 10, 'linewidth', 2);
    % Second saddle-node bifurcation of periodic orbits
plot(bifn2.snpo2.x, bifn2.snpo2.y, '-', ...
    'color', '#378c47', 'linewidth', 3);
plot(bifn2.snpo2.x, bifn2.snpo2.y, 'o', ...
    'color', '#378c47', 'markersize', 10, 'linewidth', 2);

% Format axes
xlabel("$\theta_{\mathit{u}}$", 'interpreter', 'latex');
ylabel('\tau', 'rotation', 0, 'position', [0.33, 0.25, -1], ...
    'verticalalignment', 'middle');

xlim([0.4, 1]);
xticks([0.4 0.5 0.6 0.7 0.8 0.9 1]);
xticklabels({'0.4','0.5','0.6','0.7', '0.8', '0.9', '1.0'});

ylim([0, 0.5]);
yticks([0 0.1 0.2 0.3 0.4 0.5]);
yticklabels({'0', '0.1', '0.2', '0.3', '0.4', '0.5'});

set(gca,'fontsize', 24, 'fontname', 'times');

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

p.tau_1 = 0.2; p.tau_2 = p.tau_1;
p.tspan = [0 60]; p.delays = [p.tau_1 p.tau_2];
p.options = ddeset('reltol', 1e-5);

p.history = [0.1 0.5];
[sol.dde_sim_4, nullclines] = ddeSim(p);

p.history = [0.7 0.7];
[sol.dde_sim_5] = ddeSim(p);

p.history = [0.9 0.5];
[sol.dde_sim_6] = ddeSim(p);

% Initialise figure 6
figure(6);
set(gcf, 'position', [400, 100, 600, 400]);
clf;
tiledlayout(2,2,'tilespacing','loose','padding','loose');

% Subplot 1 - \tau = 0.5
fig.ax1 = nexttile; hold on;

plot(stst.a.stable1.x, stst.a.stable1.y, '-', ...
    'color', '#c74440', 'linewidth', 3);
plot(stst.a.unstable.x, stst.a.unstable.y, '--', ...
    'color', '#c74440', 'linewidth', 3);
plot(stst.a.stable2.x, stst.a.stable2.y, '-', ...
    'color', '#c74440', 'linewidth', 3);

if isfield(po.a, 'stable1')
    plot(po.a.stable1.x, po.a.stable1.y, '-', ...
        'color', '#2c70b3', 'linewidth', 3);
    plot(po.a.stable1.x, po.a.stable1.y, 'o', ...
        'markersize', 8, 'color', '#2c70b3', 'linewidth', 2);
end
if isfield(po.a, 'unstable')
    plot(po.a.unstable.x, po.a.unstable.y, '-', ...
        'color', '#2c70b3', 'linewidth', 3);
    plot(po.a.unstable.x, po.a.unstable.y, 'x', ...
        'markersize', 8, 'color', '#2c70b3', 'linewidth', 2);
end
if isfield(po.a, 'stable2')
    plot(po.a.stable2.x, po.a.stable2.y, '-', ...
        'color', '#2c70b3', 'linewidth', 3);
    plot(po.a.stable2.x, po.a.stable2.y, 'o', ...
        'markersize', 8, 'color', '#2c70b3', 'linewidth', 2);
end

% Plot label
annotation('textbox',[0.135,0.943,0.06,0.06],'string',"(a)", ...
    'fontsize',24,'fontname','times','edgecolor','none');

% Format axes
xlabel("$\theta_{\mathit{u}}$", 'interpreter', 'latex');
ylabel("$\mathit{u}$", 'rotation', 0, 'interpreter', 'latex', ...
    'position', [0.22, 0.5, -1], 'verticalalignment', 'middle');

xlim([0.4 1]);
xticks([0.4 0.5 0.6 0.7 0.8 0.9 1]);
xticklabels({'0.4','','0.6','', '0.8', '', '1.0'});
xtickangle(0);

ylim([0 1]);
yticks([0 0.5 1]);
yticklabels({'0', '0.5', '1.0'});

set(gca,'fontsize', 24, 'fontname', 'times');

% Subplot 2 - \tau = 0.2
fig.ax2 = nexttile; hold on;

plot(stst.b.stable1.x, stst.b.stable1.y, '-', ...
    'color', '#c74440', 'linewidth', 3);
plot(stst.b.unstable.x, stst.b.unstable.y, '--', ...
    'color', '#c74440', 'linewidth', 3);
plot(stst.b.stable2.x, stst.b.stable2.y, '-', ...
    'color', '#c74440', 'linewidth', 3);

if isfield(po.b, 'stable1')
    plot(po.b.stable1.x, po.b.stable1.y, '-', ...
        'color', '#2c70b3', 'linewidth', 3);
    plot(po.b.stable1.x, po.b.stable1.y, 'o', ...
        'markersize', 8, 'color', '#2c70b3', 'linewidth', 2);
end
if isfield(po.b, 'unstable')
    plot(po.b.unstable.x, po.b.unstable.y, '-', ...
        'color', '#2c70b3', 'linewidth', 3);
    plot(po.b.unstable.x, po.b.unstable.y, 'x', ...
        'markersize', 8, 'color', '#2c70b3', 'linewidth', 2);
end
if isfield(po.b, 'stable2')
    plot(po.b.stable2.x, po.b.stable2.y, '-', ...
        'color', '#2c70b3', 'linewidth', 3);
    plot(po.b.stable2.x, po.b.stable2.y, 'o', ...
        'markersize', 8, 'color', '#2c70b3', 'linewidth', 2);
end

% Plot label
annotation('textbox',[0.593,0.943,0.06,0.06],'string',"(b)", ...
    'fontsize',24,'fontname','times','edgecolor','none');

% Format axes
xlabel("$\theta_{\mathit{u}}$", 'interpreter', 'latex');

xlim([0.5 0.9]);
xticks([0.5 0.6 0.7 0.8 0.9]);
xticklabels({'0.5','','0.7','', '0.9'});
xtickangle(0);

ylim([0 0.9]);
yticks([0.1 0.3 0.5 0.7 0.9]);
yticklabels({'0.1', '', '0.5', '', '0.9'});

set(gca,'fontsize', 24, 'fontname', 'times');

% Subplot 3 - \tau = 0.09
fig.ax3 = nexttile; hold on;

plot(stst.c.stable1.x, stst.c.stable1.y, '-', ...
    'color', '#c74440', 'linewidth', 3);
plot(stst.c.unstable.x, stst.c.unstable.y, '--', ...
    'color', '#c74440', 'linewidth', 3);
plot(stst.c.stable2.x, stst.c.stable2.y, '-', ...
    'color', '#c74440', 'linewidth', 3);

if isfield(po.c, 'stable1')
    plot(po.c.stable1.x, po.c.stable1.y, '-', ...
        'color', '#2c70b3', 'linewidth', 3);
    plot(po.c.stable1.x, po.c.stable1.y, 'o', ...
        'markersize', 8, 'color', '#2c70b3', 'linewidth', 2);
end
if isfield(po.c, 'unstable')
    plot(po.c.unstable.x, po.c.unstable.y, '-', ...
        'color', '#2c70b3', 'linewidth', 3);
    plot(po.c.unstable.x, po.c.unstable.y, 'x', ...
        'markersize', 8, 'color', '#2c70b3', 'linewidth', 2);
end
if isfield(po.c, 'stable2')
    plot(po.c.stable2.x, po.c.stable2.y, '-', ...
        'color', '#2c70b3', 'linewidth', 3);
    plot(po.c.stable2.x, po.c.stable2.y, 'o', ...
        'markersize', 8, 'color', '#2c70b3', 'linewidth', 2);
end

% Plot label
annotation('textbox',[0.135,0.465,0.06,0.06],'string',"(c)", ...
    'fontsize',24,'fontname','times','edgecolor','none');

% Format axes
xlabel("$\theta_{\mathit{u}}$", 'interpreter', 'latex');
ylabel("$\mathit{u}$", 'rotation', 0, 'interpreter', 'latex', ...
    'position', [0.384, 0.5, -1], 'verticalalignment', 'middle');

xlim([0.5 0.9]);
xticks([0.5 0.6 0.7 0.8 0.9]);
xticklabels({'0.5','','0.7', '', '0.9'});
xtickangle(0);

ylim([0.3 0.7]);
yticks([0.3 0.4 0.5 0.6 0.7]);
yticklabels({'0.3', '', '0.5', '', '0.7'});

set(gca,'fontsize', 24, 'fontname', 'times');

% Set aspect ratio 
% fig.aspectRatio = pbaspect(fig.ax1);
% pbaspect(fig.ax3, fig.aspectRatio);
% clear fig

% Subplot 4 - dde simulation
nexttile; hold on;

% Plot nullclines
plot(nullclines.u_null, nullclines.v_range, ...
    '-.', 'color', '#c74440', 'linewidth', 3);
plot(nullclines.u_range, nullclines.v_null, ...
    '-.', 'color', '#2c70b3', 'linewidth', 3);

% Plot solution
plot(sol.dde_sim_4.y(1,1), sol.dde_sim_4.y(2,1), ...
    '.', 'color', '#378c47', 'markersize', 12);
plot(sol.dde_sim_4.y(1,:), sol.dde_sim_4.y(2,:), ...
    '-', 'color', '#378c47', 'linewidth', 1.5);

plot(sol.dde_sim_5.y(1,1), sol.dde_sim_5.y(2,1), ...
    '.', 'color', '#f77e1b', 'markersize', 12);
plot(sol.dde_sim_5.y(1,:), sol.dde_sim_5.y(2,:), ...
    '-', 'color', '#f77e1b', 'linewidth', 1.5);

plot(sol.dde_sim_6.y(1,1), sol.dde_sim_6.y(2,1), ...
    '.', 'color', '#6142a6', 'markersize', 12);
plot(sol.dde_sim_6.y(1,:), sol.dde_sim_6.y(2,:), ...
    '-', 'color', '#6142a6', 'linewidth', 1.5);

% Plot label
annotation('textbox',[0.593,0.465,0.06,0.06],'string',"(d)", ...
    'fontsize',24,'fontname','times','edgecolor','none');

% Format axes
xlabel("$\mathit{u}$", 'interpreter', 'latex');
ylabel("$\mathit{v}$", 'interpreter', 'latex', 'rotation', 0, ...
    'position', [-0.275, 0.5, -1], 'verticalalignment', 'middle');

xlim([0 1]);
xticks(0:0.5:1);
xticklabels({'0','0.5','1.0'});
xtickangle(0);

ylim([0 1]);
yticks(0:0.5:1);
yticklabels({'0','0.5','1.0'});

set(gca,'fontsize', 24, 'fontname', 'times');

% Save figure
print(gcf, '../Figures/Figure_7.png', '-dpng', '-r300');

%% Calculate Maximum Lyapunov Exponents over parameter space
% This section may take up to 3 hours to complete

% Define Lyapunov exponent calculation parameters
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
p.tau_1 = 0.1; p.tau_2 = 0.1;

% Define DDE parameters
p.tspan = [0 60];
p.delays = [p.tau_1 p.tau_2];
p.history = [0.074 0.077];
p.options = ddeset('reltol', 1e-5);

% Calculate maximum Lyapunov Exponents
% le.lambda_max = calcLyapunovExponent(p);

% Alternatively, load pre-calculated lambda_max
load('lambda_max.mat');
le.lambda_max = lambda_max;
clear lambda_max

%% Figure 8 Replictae
% Extract and filter lambda_max
le.img = le.lambda_max;
le.img(le.img < -3.5) = -3.5;

% Initialise figure 7
figure(7);
set(gcf, 'position', [400, 100, 600, 400]);
clf; hold on;

% Plot lambda max
imagesc(p.a_grid(1, :), p.b_grid(:, 1), le.img);
colormap(flipud(gray));
colorbar;

% Format axes
xlabel("$\mathit{a}$", 'interpreter', 'latex');
ylabel("$\mathit{b}$", 'interpreter', 'latex', 'rotation', 0, ...
    'position', [-11, 2.5, -1], 'verticalalignment', 'middle');

xlim([-10, 0]);
xticks(-10:2:0);
xticklabels({'-10', '-8', '-6', '-4', '-2', '0'});

ylim([0, 5]);
yticks(0:1:5);
yticklabels({'0', '1', '2', '3', '4', '5'});


set(gca,'fontsize', 24, 'fontname', 'times');

% Save figure
print(gcf, '../Figures/Figure_8.png', '-dpng', '-r300');

%% Figure 9 Replicate
% Select Parameters
p.alpha = 1;
p.a = -3.7; p.b = 1.9;
p.c = 1.9; p.d = -3.7;
p.theta_u = 0.2; p.theta_v = 0.2;
p.tau_1 = 0.1; p.tau_2 = 0.1;

% Define DDE parameters
p.tspan = [0 10];
p.delays = [p.tau_1 p.tau_2];
p.options = ddeset('reltol', 1e-5);

% Define Lyapunov exponent calculation parameters
p.ptbn = 0.001; % distance to perturb points along main trajectory
p.int = 5; % interval to simulate perturbations for

% Solve DDE - chaotic solutions
p.beta = 60;
p.history = [0.1313 0.15];
sol.Chaotic_1 = ddeSim(p);
le.Chaotic_1 = lyapunovExponent(p, p.ptbn, p.int);

p.history = [0.1323 0.15];
sol.Chaotic_2 = ddeSim(p);
le.Chaotic_2 = lyapunovExponent(p, p.ptbn, p.int);

% Solve DDE - quasi-periodic solutions
p.beta = 40;
p.history = [0.1313 0.15];
sol.QuasiPeriodic_1 = ddeSim(p);
le.QuasiPeriodic_1 = lyapunovExponent(p, p.ptbn, p.int);

p.history = [0.1323 0.15];
sol.QuasiPeriodic_2 = ddeSim(p);
le.QuasiPeriodic_2 = lyapunovExponent(p, p.ptbn, p.int);

% Initialise figure 8
figure(8);
set(gcf, 'position', [400, 100, 600, 225]);
clf;
tiledlayout(1, 2, 'tilespacing', 'compact', 'padding', 'loose');

% First subplot
nexttile; hold on;

% Plot chotic solution
plot(sol.Chaotic_1.y(1,:), sol.Chaotic_1.y(2,:), ...
    '-', 'color', '#2c70b3', 'linewidth', 2);
plot(sol.Chaotic_2.y(1,:), sol.Chaotic_2.y(2,:), ...
    '-', 'color', '#f77e1b', 'linewidth', 2);

% Plot label
annotation('textbox',[0.132,0.895,0.06,0.06],'string',"(a)", ...
    'fontsize',24,'fontname','times','edgecolor','none');

% Format axes
xlabel("$\mathit{u}$", 'interpreter', 'latex');
ylabel("$\mathit{v}$", 'interpreter', 'latex', 'rotation', 0, ...
    'position', [0.078, 0.175, -1], 'verticalalignment', 'middle');

xlim([0.115, 0.235]);
xticks(0.12:0.02:0.23);
xticklabels({'', '0.14', '', '0.18', '', '0.22'});
xtickangle(0);

ylim([0.115, 0.235]);
yticks(0.12:0.02:0.23);
yticklabels({'', '0.14', '', '0.18', '', '0.22'});

set(gca,'fontsize', 24, 'fontname', 'times');

% Second subplot
nexttile; hold on;

% Plot chotic solution
plot(sol.QuasiPeriodic_1.y(1,:), sol.QuasiPeriodic_1.y(2,:), ...
    '-', 'color', '#2c70b3', 'linewidth', 2);
plot(sol.QuasiPeriodic_2.y(1,:), sol.QuasiPeriodic_2.y(2,:), ...
    '-', 'color', '#f77e1b', 'linewidth', 2);

% Plot label
annotation('textbox',[0.563,0.895,0.06,0.06],'string',"(b)", ...
    'fontsize',24,'fontname','times','edgecolor','none');

% Format axes
xlabel("$\mathit{u}$", 'interpreter', 'latex');

xlim([0.125, 0.21]);
xticks(0.13:0.01:0.21);
xticklabels({'', '0.14', '', '0.16', '', '0.18', '', '0.20', ''});
xtickangle(0);

ylim([0.125, 0.21]);
yticks(0.13:0.02:0.21);
yticklabels({'0.13', '', '0.17', '', '0.21'});

set(gca,'fontsize', 24, 'fontname', 'times');

% Save figure
print(gcf, '../Figures/Figure_9.png', '-dpng', '-r300');