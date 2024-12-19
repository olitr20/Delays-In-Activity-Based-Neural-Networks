%% DDE-BIFTOOL - Wilson-Cowan Network with Delays
%

clearvars; clc
format compact
addpath('ddebiftool/');
addpath('ddebiftool_extra_psol/');
addpath('ddebiftool_utilities/');
%#ok<*ASGLU,*NOPTS,*NASGU>

%% Definition of User Functions
% Right-hand side
neuron_sys_rhs=@(xx,par)[...
    -xx(1,1)+1/(1+exp(-par(2)*(par(7)+par(3)*xx(1,2)+par(4)*xx(2,2))));....
    par(1)*(-xx(2,1)+1/(1+exp(-par(2)*(par(8)+par(5)*xx(1,2)+par(6)*xx(2,2)))))];

% Delays and Continuation Parameters
neuron_tau = @() 9;
ind_theta_u = 7;
ind_taus = 9;

funcs=set_funcs(...
    'sys_rhs', neuron_sys_rhs,...
    'sys_tau', @() 9);

% Define model parameters
p.alpha = 1; p.beta = 60;
p.a = -1; p.b = -0.4; p.c = -1; p.d = 0;
p.theta_u = 0.7; p.theta_v = 0.5;
p.tau_1 = 0.09; p.tau_2 = p.tau_1;

% Define continuation parameters
if p.tau_1 == 0.09
    min_theta_u = 0.61;
    max_theta_u = 1;
elseif p.tau_1 == 0.2
    min_theta_u = 0.4;
    max_theta_u = 0.85;
elseif p.tau_1 == 0.5
    min_theta_u = 0.4;
    max_theta_u = 0.97;
else
    min_theta_u = 0.4;
    max_theta_u = 1;
end

%% Initial Guess for Steady State
stst.kind='stst';
stst.parameter=[p.alpha, p.beta, p.a, p.b, p.c, p.d, ...
    p.theta_u, p.theta_v, p.tau_1];

% stst.x=[0.103596478880629;1]; % theta_u = 0.5
% stst.x=[0.599328884047436;0]; % theta_u = 0.6
stst.x=[0.698598941633828;0]; % theta_u = 0.7
% stst.x=[0.797713228634354;0]; % theta_u = 0.8
% stst.x=[0.896403521119372;0]; % theta_u = 0.9

method = df_mthod(funcs, 'stst'); % flag_newhheur omitted
method.stability.minimal_real_part = -2;

% Correct steady state point
[stst, success] = p_correc(funcs, stst, [], [], method.point);
if success == 1
    disp('Initial Steady State Correction Successful')
end

% Compute its stability
stst.stability = p_stabil(funcs, stst, method.stability);

%% Initialize branch of steady state
% Get an empty branch with theta_u as a free parameter
branch8 = df_brnch(funcs, ind_theta_u, 'stst');

% Set bounds for continuation parameter
branch8.parameter.min_bound(1,:) = [ind_theta_u 0.4];
branch8.parameter.max_bound(1,:) = [ind_theta_u 1];
branch8.parameter.max_step(1,:) = [ind_theta_u 0.005];

% Use steady state point as first steady state branch point
branch8.point = stst;

%% Continue branch of steady state in theta_u
% Perturb steady state point by 0.01 in theta_u
stst.parameter(ind_theta_u) = stst.parameter(ind_theta_u) + 0.005;

% Correct new steady state point
[stst, success] = p_correc(funcs, stst, [], [], method.point);
if success == 1
    disp('Perturbed Steady State Correction Successful')
end

% Use new steady state point as second steady state branch point
branch8.point(2) = stst;
branch8.method.continuation.plot = 1;

% Initialise figure 1
figure(1); clf;

% Continue steady state branch in both directions
branch8 = br_contn(funcs, branch8, 200);
branch8 = br_rvers(branch8);
branch8 = br_contn(funcs, branch8, 200);

% Format figure 1
xlabel("$\theta_{\mathit{u}}$", 'Interpreter', 'latex');
ylabel("$\mathit{u}$", 'Rotation', 0, 'Interpreter', 'latex');
xlim([0.4 1]);
ylim([0 1]);

%% Compute Stability of this Steady State Branch
% Calculate stability of every point along branch
branch8.method.stability.minimal_real_part = -2;
branch8 = br_stabl(funcs, branch8, 0, 0);

% Obtain suitable scalar measures to plot stability along branch
[xm, ym] = df_measr(1, branch8);
ym.subfield = 'l0';

% Initialise figure 2
figure(2); clf;

% Plot stability along branch
br_plot(branch8, [], ym, 'b');
br_plot(branch8, [], ym, 'b.');
plot([0 length(branch8.point)], [0 0], '-.');

% Format figure 2
xlabel('Point Number Along Branch');
ylabel('\Re(\lambda)', 'Rotation', 0);
xlim([0 length(branch8.point)]);

%% Locate the First Hopf Point
% Find the hopf bifurcation point
ind_hopf = find(arrayfun(@(x) real(x.stability.l0(1)) > 0, ...
    branch8.point), 1, 'last');

% Convert this point to a hopf point
hopf = p_tohopf(funcs, branch8.point(ind_hopf));

% Get hopf calculation method parameters
method = df_mthod(funcs, 'hopf');
method.stability.minimal_real_part = -1;

% Correct hopf point
[hopf, success] = p_correc(funcs, hopf, ind_theta_u, [], method.point);
if success == 1
    disp('Initial Hopf Point Correction Successful')
end
first_hopf = hopf; % store hopf point for later use

%% Initialise and Continue Periodic Orbits
% Using the first Hopf point, construct a small-amplitude (1e-2) periodic
% solution on an equidistant mesh of 18 intervals with piecewise polynomial
% degree 3. The steplength condition (returned by p_topsol) ensures the
% branch switch from the Hopf to the periodic solution as it avoids
% convergence of the amplitude to zero during corrections.

% Construct small amplitude oscillation
intervals = 18;
degree = 3;
[psol, stepcond] = p_topsol(funcs, first_hopf, 1e-2, degree, intervals);

% Correct periodic solution
method = df_mthod(funcs, 'psol');
[psol, success] = p_correc(funcs, psol, ind_theta_u, stepcond, method.point);
if success == 1
    disp('Periodic Solution Correction Successful')
end

% Get an empty branch with theta_u as a free parameter
branch9=df_brnch(funcs,ind_theta_u,'psol');

% Set bounds for continuation parameter
branch9.parameter.min_bound(1,:) = [ind_theta_u min_theta_u];
branch9.parameter.max_bound(1,:) = [ind_theta_u max_theta_u];
branch9.parameter.max_step(1,:) = [ind_theta_u 0.005];

% Make degenerate periodic solution with amplitude zero at hopf point
deg_psol = p_topsol(funcs, first_hopf, 0, degree, intervals);

% Use deg_psol and psol as first two points on branch
deg_psol.mesh = []; % clear the mesh field to save memory and avoid adaptive mesh selection
branch9.point = deg_psol;
psol.mesh = [];
branch9.point(2) = psol;

% Initialise figure 3
figure(3); clf;

% Continue periodic solutions branch, plotting amplitude
branch9 = br_contn(funcs, branch9, 200);

% Format figure 3
xlabel("$\theta_{\mathit{u}}$", 'Interpreter', 'latex');
ylabel('Amplitude');

%% Extract max u and theta_u along the branch 
[xm, ym] = df_measr(0, branch9);
ym.field = 'profile'; ym.row = 1; ym.col = 'all';
po_y = br_measr(branch9,ym);
po_y = max(po_y, [], 2);
po_x = br_measr(branch9,xm);

% Initialise figure 4
figure(1);
hold on

plot(po_x, po_y, 'r-')
plot(po_x, po_y, 'r.')

%% Compute Stability of this Periodic Orbit Branch
% Calculate stability of every point along branch
branch9.method.stability.minimal_real_part = -2;
branch9 = br_stabl(funcs, branch9, 0, 0);

% Obtain suitable scalar measures to plot stability along branch
[xm, ym] = df_measr(1, branch9);
ym.subfield = 'mu'; ym.row = 1;

% Initialise figure 4
figure(4); clf;

% Plot stability along branch
br_plot(branch9, xm, ym, 'b');
br_plot(branch9, xm, ym, 'b.');
% plot([0 length(branch8.point)], [0 0], '-.');

% Format figure 4
xlabel('Point Number Along Branch');
ylabel('\Re(\lambda)', 'Rotation', 0);

%% Construct Final Figure
% Initialise figure 5
figure(5); clf
hold on

[xm, ym] = df_measr(0, branch8);
stst_x = br_measr(branch8, xm);
stst_y = br_measr(branch8, ym);

ind_hopf_1 = find(arrayfun(@(x) real(x.stability.l0(1)) > 0, ...
    branch8.point), 1, 'first');
ind_hopf_2 = find(arrayfun(@(x) real(x.stability.l0(1)) > 0, ...
    branch8.point), 1, 'last');
if isempty(ind_hopf_1)
    ind_hopf_1 = length(branch8.point);
end

stst_x_stable_1 = stst_x(1:ind_hopf_1);
stst_y_stable_1 = stst_y(1:ind_hopf_1);
stst_x_unstable = stst_x(ind_hopf_1:ind_hopf_2);
stst_y_unstable = stst_y(ind_hopf_1:ind_hopf_2);
stst_x_stable_2 = stst_x(ind_hopf_2:end);
stst_y_stable_2 = stst_y(ind_hopf_2:end);

plot(stst_x_stable_1, stst_y_stable_1, 'k-')
plot(stst_x_unstable, stst_y_unstable, 'k--')
plot(stst_x_stable_2, stst_y_stable_2, 'k-')

[xm, ym] = df_measr(0, branch9);
ym.field = 'profile'; ym.row = 1; ym.col = 'all';
po_x = br_measr(branch9,xm);
po_y = br_measr(branch9,ym);
po_y = max(po_y, [], 2);

ind_hopf_po_1 = find(arrayfun(@(x) real(x.stability.mu(1)) > 1.1, ...
    branch9.point), 1, 'first');
ind_hopf_po_2 = find(arrayfun(@(x) real(x.stability.mu(1)) > 1.1, ...
    branch9.point), 1, 'last');
if isempty(ind_hopf_po_1)
    ind_hopf_po_1 = length(branch9.point);
end


po_x_stable_1 = po_x(1:ind_hopf_po_1);
po_y_stable_1 = po_y(1:ind_hopf_po_1);
po_x_unstable = po_x(ind_hopf_po_1:ind_hopf_po_2);
po_y_unstable = po_y(ind_hopf_po_1:ind_hopf_po_2);
po_x_stable_2 = po_x(ind_hopf_po_2:end);
po_y_stable_2 = po_y(ind_hopf_po_2:end);


if ind_hopf_po_1 ~= 1
    plot(po_x_stable_1, po_y_stable_1, 'k-o', 'markersize', 4)
end

plot(po_x_unstable, po_y_unstable, 'k-x', 'markersize', 4)

if ind_hopf_po_2 ~= length(branch9.point)
    plot(po_x_stable_2, po_y_stable_2, 'k-o', 'markersize', 4)
end

% Format figure 5
xlabel("$\theta_{\mathit{u}}$", 'Interpreter', 'latex');
ylabel("$\mathit{u}$", 'Rotation', 0, 'Interpreter', 'latex');
xlim([0.4 1]);
ylim([0 1]);
xticks([0.4 0.5 0.6 0.7 0.8 0.9 1])
xticklabels({'0.4','0.5','0.6','0.7', '0.8', '0.9', '1.0'})
yticks([0 0.2 0.4 0.6 0.8 1.0])
yticklabels({'0', '0.2', '0.4', '0.6', '0.8', '1.0'});

% stop br_contn when matrix close to singular
% sample branches at regular theta_us to even out