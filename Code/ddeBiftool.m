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
method = df_mthod(funcs, 'hopf'); % flag_newhheur omitted
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
branch2 = df_brnch(funcs, [ind_theta_u, ind_taus], 'hopf');

% Set bounds for continuation parameters
branch2.parameter.min_bound(1:2,:) = [[ind_theta_u 0.4]' [ind_taus 0]']';
branch2.parameter.max_bound(1:2,:) = [[ind_theta_u 1]' [ind_taus 0.5]']';
branch2.parameter.max_step(1:2,:) = [[ind_theta_u 0.005]' [ind_taus 0.005]']';

% Use hopf point as first hopf branch point
branch2.point = hopf;

% Perturb hopf point by 0.0001 in theta_u
hopf.parameter(ind_theta_u) = hopf.parameter(ind_theta_u) - 0.0001;

% Correct hopf point in tau and recompute stability
[hopf, success] = p_correc(funcs, hopf, ind_taus, [], method.point);
if success == 1
    disp('Perturbed Hopf Point Correction Successful')
end

% Use perturbed hopf point as second hopf branch point
branch2.point(2) = hopf;

% Initialise figure 4
figure(4); clf;

% Continue hopf branch in both directions
branch2 = br_contn(funcs, branch2, 500);
branch2 = br_rvers(branch2);
branch2 = br_contn(funcs, branch2, 500);

% Format figure 4
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
branch3=df_brnch(funcs,ind_theta_u,'psol');

% Set bounds for continuation parameter
branch3.parameter.min_bound(1,:) = [ind_theta_u 0];
branch3.parameter.max_bound(1,:) = [ind_theta_u 0.84];
branch3.parameter.max_step(1,:) = [ind_theta_u 0.005];

% Make degenerate periodic solution with amplitude zero at hopf point
deg_psol = p_topsol(funcs, first_hopf, 0, degree, intervals);

% Use deg_psol and psol as first two points on branch
deg_psol.mesh = []; % clear the mesh field to save memory and avoid adaptive mesh selection
branch3.point = deg_psol;
psol.mesh = [];
branch3.point(2) = psol;

% Initialise figure 5
figure(5); clf;

% Continue periodic solutions branch, plotting amplitude
branch3 = br_contn(funcs, branch3, 500);

% Format figure 5
xlabel("$\theta_{\mathit{u}}$", 'Interpreter', 'latex');
ylabel('Amplitude');

% Initialise figure 6
figure(6); clf;

% Plot amplitude along the branch
[xm, ym] = df_measr(0, branch3); 
ym.field = 'period'; ym.col = 1;
br_plot(branch3,xm,ym,'b');
axis([0.545 0.855 0.735 0.775]);

% Format figure 6
xlabel("$\theta_{\mathit{u}}$", 'Interpreter', 'latex');
ylabel('Period');

%% Initialise Folds of Periodic Orbits
% something going wrong here %

% Speed up computations by vectorisation
neuron_sys_rhs = @(xx,par) [...
    -xx(1,1,:)+1/(1+exp(-par(2)*(par(7)+par(3)*xx(1,2,:)+par(4)*xx(2,2,:))));....
    par(1)*(-xx(2,1,:)+1/(1+exp(-par(2)*(par(8)+par(5)*xx(1,2,:)+par(6)*xx(2,2,:)))))];
vfuncs=set_funcs(...
    'sys_rhs', neuron_sys_rhs,...
    'sys_tau', @() 9,...
    'x_vectorized', true);

% Find initial guess
[~, indmax1] = max(arrayfun(@(x) x.parameter(ind_theta_u), branch3.point));
[~, indmax2] = min(arrayfun(@(x) x.parameter(ind_theta_u), branch3.point));

% Initialise branch and set up extended system
branch3.method.point.newton_max_iterations = 16;
[foldfuncs, branch4] = SetupPOfold(vfuncs, branch3, indmax1, ...
    'contpar', [ind_theta_u, ind_taus], 'dir', ind_taus, ...
    'print_residual_info', 1, 'step', 0.01, 'plot_measure', [],...
    'min_bound', [ind_theta_u,0.572535; ind_taus,0.102212], ...
    'max_bound', [ind_theta_u,0.827663; ind_taus,0.5],...
    'max_step', [ind_theta_u,0.01; ind_taus,0.01]);

%% Continue Folds of Periodic Orbits
% Initialise figure 7
figure(7); clf;

% Continue one fold
branch4.method.point.print_residual_info = 0;
branch4 = br_contn(foldfuncs, branch4, 100);
branch4 = br_rvers(branch4);
branch4 = br_contn(foldfuncs, branch4, 100);

% Initilaise and continue second periodic orbit branch
[foldfuncs, branch5] = SetupPOfold(vfuncs, branch3, indmax2, ...
    'contpar', [ind_theta_u, ind_taus], 'dir', ind_taus, ...
    'print_residual_info', 1, 'step', 0.01, 'plot_measure', [],...
    'min_bound', [ind_theta_u,0.4; ind_taus,0.0], ...
    'max_bound', [ind_theta_u,1.0; ind_taus,0.5],...
    'max_step', [ind_theta_u,0.01; ind_taus,0.01]);

% Continue second fold
branch5.method.point.print_residual_info = 0;
branch5 = br_contn(foldfuncs, branch5, 100);
branch5 = br_rvers(branch5);
branch5 = br_contn(foldfuncs, branch5, 100);

% Format figure 7
xlabel("$\theta_{\mathit{u}}$", 'Interpreter', 'latex');
ylabel('\tau');