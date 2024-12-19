%% DDE-BIFTOOL - Wilson-Cowan Network with Delays
% Could convert this into a function to be called from delay_wilson_cowan.m

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
p.tau_1 = 0.2; p.tau_2 = p.tau_1;

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
branch1 = df_brnch(funcs, ind_theta_u, 'stst');

% Set bounds for continuation parameter
branch1.parameter.min_bound(1,:) = [ind_theta_u 0.4];
branch1.parameter.max_bound(1,:) = [ind_theta_u 1];
branch1.parameter.max_step(1,:) = [ind_theta_u 0.005];

% Use steady state point as first steady state branch point
branch1.point = stst;

%% Continue branch of steady state in theta_u
% Perturb steady state point by 0.01 in theta_u
stst.parameter(ind_theta_u) = stst.parameter(ind_theta_u) + 0.005;

% Correct new steady state point
[stst, success] = p_correc(funcs, stst, [], [], method.point);
if success == 1
    disp('Perturbed Steady State Correction Successful')
end

% Use new steady state point as second steady state branch point
branch1.point(2) = stst;
branch1.method.continuation.plot = 1;

% Initialise figure 1
figure(1); clf;

% Continue steady state branch in both directions
branch1 = br_contn(funcs, branch1, 200);
branch1 = br_rvers(branch1);
branch1 = br_contn(funcs, branch1, 200);

% Format figure 1
xlabel("$\theta_{\mathit{u}}$", 'Interpreter', 'latex');
ylabel("$\mathit{u}$", 'Rotation', 0, 'Interpreter', 'latex');

%% Compute Stability of this Steady State Branch
% Calculate stability of every point along branch
branch1.method.stability.minimal_real_part = -2;
branch1 = br_stabl(funcs, branch1, 0, 0);

% Obtain suitable scalar measures to plot stability along branch
[xm, ym] = df_measr(1, branch1);
ym.subfield = 'l0';

% Initialise figure 2
figure(2); clf;

% Plot stability along branch
br_plot(branch1, [], ym, 'b');
br_plot(branch1, [], ym, 'b.');
plot([0 length(branch1.point)], [0 0], '-.');

% Format figure 2
xlabel('Point Number Along Branch');
ylabel('\Re(\lambda)', 'Rotation', 0);
xlim([0 length(branch1.point)]);

%% Locate the First Saddle Node Point
% branch1.point(88) and branch1.point(141)
x1s = arrayfun(@(x) x.x(1), branch1.point);
x2s = arrayfun(@(x) x.x(2), branch1.point);

figure(3); clf; hold on
plot(x1s)
plot(x2s)

% saddle nodes occur when there is one 0 eigenvalue?

% Find the first saddle-node bifurcation point
ind_sn = 88;

% Convert this point to a saddle-node point
sn = p_tofold(funcs, branch1.point(ind_sn));

% Get saddle-node calculation method parameters
method = df_mthod(funcs, 'fold');
method.stability.minimal_real_part = -1;

% Correct saddle-node point
[sn, success] = p_correc(funcs, sn, ind_theta_u, [], method.point);
if success == 1
    disp('Initial First Saddle-Node Point Correction Successful')
end

% Compute stability of saddle-node point
sn.stability = p_stabil(funcs, sn, method.stability);

%% Initialise and Continue Saddle-Node Bifurcation
% Get an empty branch with theta_u and tau as free parameters
branch2 = df_brnch(funcs, [ind_theta_u, ind_taus], 'fold');

% Set bounds for continuation parameters
branch2.parameter.min_bound(1:2,:) = [[ind_theta_u 0.4]' [ind_taus 0]']';
branch2.parameter.max_bound(1:2,:) = [[ind_theta_u 1]' [ind_taus 0.5]']';
branch2.parameter.max_step(1:2,:) = [[ind_theta_u 0.005]' [ind_taus 0.005]']';

% Use saddle-node point as first saddle-node branch point
branch2.point = sn;

% Perturb saddle-node point by 0.0001 in tau
sn.parameter(ind_taus) = sn.parameter(ind_taus) + 0.0001;

% Correct saddle-node point in theta_u and recompute stability
[sn, success] = p_correc(funcs, sn, ind_theta_u, [], method.point);
if success == 1
    disp('Perturbed First Saddle-Node Point Correction Successful')
end

% Use perturbed hopf point as second hopf branch point
branch2.point(2) = sn;

% Initialise figure 4
figure(4); clf;

% Continue hopf branch in both directions
branch2 = br_contn(funcs, branch2, 500);
branch2 = br_rvers(branch2);
branch2 = br_contn(funcs, branch2, 500);

% Format figure 4
xlabel("$\theta_{\mathit{u}}$", 'Interpreter', 'latex');
ylabel('\tau', 'Rotation', 0);

%% Find and Continue Second Saddle-Node Bifurcation
% Find the second saddle-node bifurcation point
ind_sn = 141;

% Convert this point to a saddle-node point
sn = p_tofold(funcs, branch1.point(ind_sn));

% Correct saddle-node point
[sn, success] = p_correc(funcs, sn, ind_theta_u, [], method.point);
if success == 1
    disp('Initial Second Saddle-Node Point Correction Successful')
end

% Compute stability of saddle-node point
sn.stability = p_stabil(funcs, sn, method.stability);

% Get an empty branch with theta_u and tau as free parameters
branch3 = df_brnch(funcs, [ind_theta_u, ind_taus], 'fold');

% Set bounds for continuation parameters
branch3.parameter.min_bound(1:2,:) = [[ind_theta_u 0.4]' [ind_taus 0]']';
branch3.parameter.max_bound(1:2,:) = [[ind_theta_u 1]' [ind_taus 0.5]']';
branch3.parameter.max_step(1:2,:) = [[ind_theta_u 0.005]' [ind_taus 0.005]']';

% Use saddle-node point as first saddle-node branch point
branch3.point = sn;

% Perturb saddle-node point by 0.0001 in tau
sn.parameter(ind_taus) = sn.parameter(ind_taus) + 0.0001;

% Correct saddle-node point in theta_u and recompute stability
[sn, success] = p_correc(funcs, sn, ind_theta_u, [], method.point);
if success == 1
    disp('Perturbed Second Saddle-Node Point Correction Successful')
end

% Use perturbed hopf point as second hopf branch point
branch3.point(2) = sn;

% Continue hopf branch in both directions
branch3 = br_contn(funcs, branch3, 500);
branch3 = br_rvers(branch3);
branch3 = br_contn(funcs, branch3, 500);

% Format figure 4
xlabel("$\theta_{\mathit{u}}$", 'Interpreter', 'latex');
ylabel('\tau', 'Rotation', 0);
ylim([0 0.5])

%% Locate the First Hopf Point
% Where eigenvalue curves in the stability plot (figure 2) cross the zero 
% line, hopf bifurcations occur. To find this hopf point, select the last
% point with positive eigenvalues and turn it into an (approximate) Hopf
% bifurcation point. Then, correct the Hopf point using appropriate method
% parameters and one free parameter (theta_u).

% Find the hopf bifurcation point
ind_hopf = find(arrayfun(@(x) real(x.stability.l0(1)) > 0, ...
    branch1.point), 1, 'last');

% Convert this point to a hopf point
hopf = p_tohopf(funcs, branch1.point(ind_hopf));

% Get hopf calculation method parameters
method = df_mthod(funcs, 'hopf');
method.stability.minimal_real_part = -1;

% Correct hopf point
[hopf, success] = p_correc(funcs, hopf, ind_theta_u, [], method.point);
if success == 1
    disp('Initial Hopf Point Correction Successful')
end
first_hopf = hopf; % store hopf point for later use

% Compute stability of hopf point
hopf.stability = p_stabil(funcs, hopf, method.stability);

%% Initialise and Continue Hopf Bifurcation
% In order to follow a branch of Hopf bifurcations in the two parameter
% space, (theta_u, tau), two starting points are, again, required. Hence
% the Hopf point already found is used as well as one perturbed in tau and 
% corrected in theta_u. For the free parameters, theta_u and tau, suitable
% intervals and maximal step sizes are provided.

% Get an empty branch with theta_u and tau as free parameters
branch4 = df_brnch(funcs, [ind_theta_u, ind_taus], 'hopf');

% Set bounds for continuation parameters
branch4.parameter.min_bound(1:2,:) = [[ind_theta_u 0.4]' [ind_taus 0]']';
branch4.parameter.max_bound(1:2,:) = [[ind_theta_u 1]' [ind_taus 0.5]']';
branch4.parameter.max_step(1:2,:) = [[ind_theta_u 0.005]' [ind_taus 0.005]']';

% Use hopf point as first hopf branch point
branch4.point = hopf;

% Perturb hopf point by 0.0001 in theta_u
hopf.parameter(ind_theta_u) = hopf.parameter(ind_theta_u) - 0.0001;

% Correct hopf point in tau and recompute stability
[hopf, success] = p_correc(funcs, hopf, ind_taus, [], method.point);
if success == 1
    disp('Perturbed Hopf Point Correction Successful')
end

% Use perturbed hopf point as second hopf branch point
branch4.point(2) = hopf;

% Initialise figure 5
figure(5); clf;

% Continue hopf branch in both directions
branch4 = br_contn(funcs, branch4, 500);
branch4 = br_rvers(branch4);
branch4 = br_contn(funcs, branch4, 500);

% Format figure 5
xlabel("$\theta_{\mathit{u}}$", 'Interpreter', 'latex');
ylabel('\tau', 'Rotation', 0);

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
branch5=df_brnch(funcs,ind_theta_u,'psol');

% Set bounds for continuation parameter
branch5.parameter.min_bound(1,:) = [ind_theta_u 0];
branch5.parameter.max_bound(1,:) = [ind_theta_u 0.85];
branch5.parameter.max_step(1,:) = [ind_theta_u 0.005];

% Make degenerate periodic solution with amplitude zero at hopf point
deg_psol = p_topsol(funcs, first_hopf, 0, degree, intervals);

% Use deg_psol and psol as first two points on branch
deg_psol.mesh = []; % clear the mesh field to save memory and avoid adaptive mesh selection
branch5.point = deg_psol;
psol.mesh = [];
branch5.point(2) = psol;

% Initialise figure 6
figure(6); clf;

% Continue periodic solutions branch, plotting amplitude
branch5 = br_contn(funcs, branch5, 500);

% Format figure 6
xlabel("$\theta_{\mathit{u}}$", 'Interpreter', 'latex');
ylabel('Amplitude');
xlim([0.55 0.85]);

%% Initialise Folds of Periodic Orbits
% Speed up computations by vectorisation
neuron_sys_rhs = @(xx,par) [...
    -xx(1,1,:)+1/(1+exp(-par(2)*(par(7)+par(3)*xx(1,2,:)+par(4)*xx(2,2,:))));....
    par(1)*(-xx(2,1,:)+1/(1+exp(-par(2)*(par(8)+par(5)*xx(1,2,:)+par(6)*xx(2,2,:)))))];
vfuncs=set_funcs(...
    'sys_rhs', neuron_sys_rhs,...
    'sys_tau', @() 9,...
    'x_vectorized', true);

% Extract amplitude along the branch 
[~, ym] = df_measr(0, branch5);
ym.field = 'profile'; ym.row = 1; ym.col = 'ampl';
amplitudes = br_measr(branch5,ym);
val_ind = find(amplitudes > 0.1); % restrict branch to ampl > 0.1

% Find the maximum theta_u where amplitudes > 0.1
[~, local_indmax] = max(arrayfun(@(x) x.parameter(ind_theta_u), ...
    branch5.point(val_ind)));
indmax = val_ind(local_indmax);

% Find the maximum theta_u where amplitudes > 0.1
[~, local_indmin] = min(arrayfun(@(x) x.parameter(ind_theta_u), ...
    branch5.point(val_ind)));
indmin = val_ind(local_indmin);

%% Initialise and Continue First Fold of Periodic Orbits
psol = branch5.point(indmax:indmax+1);
intervals = 40; degree = 4;
psol = arrayfun(@(p) p_remesh(p, degree, intervals), psol); % refine
method.point.adapt_mesh_after_correct = 1;
method.point.newton_max_iterations = 7;
method.point.newton_nmon_iterations = 2;
psol = arrayfun(@(p) p_correc(funcs, p, [], [], method.point), psol); %correct

% Copy branch 5 and start at max theta_u point
branch6 = branch5;
branch6.point = psol;
branch6.method = method;

% Initialise and set up branch
branch6.method.point.newton_max_iterations = 16;
[foldfuncs, branch6] = SetupPOfold(vfuncs, branch6, 1, ...
    'contpar', [ind_theta_u, ind_taus], 'dir', ind_taus, ...
    'print_residual_info', 1, 'step', 0.01, 'plot_measure', [],...
    'min_bound', [ind_theta_u,0.572536; ind_taus,0.102213], ...
    'max_bound', [ind_theta_u,0.827662; ind_taus,0.5],...
    'max_step', [ind_theta_u,0.01; ind_taus,0.01]);

% Initialise figure 7
figure(7); clf;

% Continue one fold
branch6.method.point.print_residual_info = 0;
branch6 = br_contn(foldfuncs, branch6, 100);
branch6 = br_rvers(branch6);
branch6 = br_contn(foldfuncs, branch6, 100);

%% Initialise and Continue Second Fold of Periodic Orbits
psol = branch5.point(indmin-1:indmin);
intervals = 40; degree = 4;
psol = arrayfun(@(p) p_remesh(p, degree, intervals), psol); % refine
method.point.adapt_mesh_after_correct = 1;
method.point.newton_max_iterations = 7;
method.point.newton_nmon_iterations = 2;
psol = arrayfun(@(p) p_correc(funcs, p, [], [], method.point), psol); %correct

% Copy branch 5 and start at max theta_u point
branch7 = branch5;
branch7.point = psol;
branch7.method = method;

% Initilaise and continue second periodic orbit branch
[foldfuncs, branch7] = SetupPOfold(vfuncs, branch7, 1, ...
    'contpar', [ind_theta_u, ind_taus], 'dir', ind_taus, ...
    'print_residual_info', 1, 'step', 0.01, 'plot_measure', [],...
    'min_bound', [ind_theta_u,0.4; ind_taus,0.0], ...
    'max_bound', [ind_theta_u,1.0; ind_taus,0.5],...
    'max_step', [ind_theta_u,0.01; ind_taus,0.01]);

% Continue second fold
branch7.method.point.print_residual_info = 0;
branch7 = br_contn(foldfuncs, branch7, 100);
branch7 = br_rvers(branch7);
branch7 = br_contn(foldfuncs, branch7, 100);

% Format figure 7
xlabel("$\theta_{\mathit{u}}$", 'Interpreter', 'latex');
ylabel('\tau');

%% Construct Final Figure
% Get deafult plotting measures: x = theta_u, y = tau
[xm, ym] = df_measr(0, branch2);

% Extract data from branches
% % potentially use deval to evaluate each branch at a grid of points to
% % even out the lines.

    % First saddle-node bifurcation
    sn_bifn_1_x = br_measr(branch2, xm);
    sn_bifn_1_y = br_measr(branch2, ym);

    % Second saddle-node bifurcation
    sn_bifn_2_x = br_measr(branch3, xm);
    sn_bifn_2_y = br_measr(branch3, ym);

    % Hopf bifurcation
    hopf_bifn_x = br_measr(branch4, xm);
    hopf_bifn_y = br_measr(branch4, ym);

    % First saddle-node of periodic orbits
    snpo_bifn_1_x = br_measr(branch6, xm);
    snpo_bifn_1_x = linspace(max(snpo_bifn_1_x), min(snpo_bifn_1_x), 19);
    snpo_bifn_1_x = snpo_bifn_1_x(1:12);
    snpo_bifn_1_y = br_measr(branch6, ym);
    snpo_bifn_1_y = linspace(min(snpo_bifn_1_y), max(snpo_bifn_1_y), 19);
    snpo_bifn_1_y = snpo_bifn_1_y(1:12);

    % Second saddle-node of periodic orbits
    snpo_bifn_2_x = br_measr(branch7, xm);
    snpo_bifn_2_x = linspace(min(snpo_bifn_2_x), max(snpo_bifn_2_x), 19);
    snpo_bifn_2_x = snpo_bifn_2_x(1:12);
    snpo_bifn_2_y = br_measr(branch7, ym);
    snpo_bifn_2_y = linspace(min(snpo_bifn_2_y), max(snpo_bifn_2_y), 19);
    snpo_bifn_2_y = snpo_bifn_2_y(1:12);

% Initialise figure 8
figure(8); clf;
hold on

% Plot bifurcations
plot(sn_bifn_1_x, sn_bifn_1_y, 'k', 'linewidth', 1)
plot(sn_bifn_2_x, sn_bifn_2_y, 'k', 'linewidth', 1)
dashline(hopf_bifn_x, hopf_bifn_y, ...
    1.5, 1, 1.5, 1, 'color', 'k', 'linewidth', 1)
plot(snpo_bifn_1_x, snpo_bifn_1_y, 'k-o', 'markersize', 4, 'linewidth', 1)
plot(snpo_bifn_2_x, snpo_bifn_2_y, 'k-o', 'markersize', 4, 'linewidth', 1)

% Format figure 8
set(gca,'FontSize', 14, 'FontName', 'Times')
xlabel("$\theta_{\mathit{u}}$", 'Interpreter', 'latex')
ylabel('\tau', 'rotation', 0);
xlim([0.4, 1]);
ylim([0, 0.5]);
xticks([0.4 0.5 0.6 0.7 0.8 0.9 1])
xticklabels({'0.4','0.5','0.6','0.7', '0.8', '0.9', '1.0'})
yticks([0 0.1 0.2 0.3 0.4 0.5])
yticklabels({'0', '0.1', '0.2', '0.3', '0.5', '0.5'});

