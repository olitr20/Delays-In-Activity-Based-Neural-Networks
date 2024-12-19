%% DDE-BIFTOOL - Wilson-Cowan Network with Delays

function bifn = ddeBiftoolMain(p)

addpath('ddebiftool/');
addpath('ddebiftool_extra_psol/');
addpath('ddebiftool_utilities/');

%#ok<*ASGLU,*NOPTS,*NASGU>

%% Define Model
% Right-hand side
sys_rhs = @(xx, par) [ ...
        -xx(1,1) + 1 / (1 + exp(-par(2) * (par(7) + par(3) * xx(1,2) + par(4) * xx(2,2)))); ...
        par(1) * (-xx(2,1) + 1 / (1 + exp(-par(2) * (par(8) + par(5) * xx(1,2) + par(6) * xx(2,2))))) ...
    ];

% Delays and Continuation Parameters
sys_tau = @() 9;
ind_theta_u = 7;
ind_taus = 9;

funcs=set_funcs( ...
    'sys_rhs', sys_rhs, ...
    'sys_tau', sys_tau);

%% Initial Guess for Steady State
% Initialise steady state point
stst.kind='stst';
stst.parameter=[p.alpha, p.beta, p.a, p.b, p.c, p.d, ...
    p.theta_u, p.theta_v, p.tau_1];
stst.x=[0.698598941633828;0]; % Extract [u, v] from odeSim endpoint

% Initialise point correction method
method = df_mthod(funcs, 'stst');
method.stability.minimal_real_part = -2;

% Correct steady state point
[stst, success] = p_correc(funcs, stst, [], [], method.point);
if success ~= 1
    disp('Initial Steady State Correction Unsuccessful')
end

% Compute point stability
stst.stability = p_stabil(funcs, stst, method.stability);

%% Initialise steady state branch
% Get an empty branch with theta_u as a free parameter
branch1 = df_brnch(funcs, ind_theta_u, 'stst');

% Set bounds for continuation parameter
branch1.parameter.min_bound(1,:) = [ind_theta_u 0.4];
branch1.parameter.max_bound(1,:) = [ind_theta_u 1];
branch1.parameter.max_step(1,:) = [ind_theta_u 0.005];

% Use steady state point as first branch point
branch1.point = stst;

%% Continue branch of steady state in theta_u
% Perturb steady state point in theta_u
stst.parameter(ind_theta_u) = stst.parameter(ind_theta_u) + 0.005;

% Correct new steady state point
[stst, success] = p_correc(funcs, stst, [], [], method.point);
if success ~= 1
    disp('Perturbed Steady State Correction Unsuccessful')
end

% Use new steady state point as second branch point
branch1.point(2) = stst;

% Continue steady state branch in both directions
branch1.method.continuation.plot = 0; % hide continuation plot
branch1 = br_contn(funcs, branch1, 200);
branch1 = br_rvers(branch1);
branch1 = br_contn(funcs, branch1, 200);

% Calculate stability along branch
branch1.method.stability.minimal_real_part = -2;
branch1 = br_stabl(funcs, branch1, 0, 0);

%% Locate First Saddle Node Point
[xm, ym] = df_measr(0, branch1); % extract default measures for branch 1
us = br_measr(branch1, ym); % extract u values from branch 1

% Find the minimum theta_u where us > 0.5
val_ind = find(us > 0.5); % restrict branch to u > 0.5
[~, local_indmin] = min(arrayfun(@(x) x.parameter(ind_theta_u), ...
    branch1.point(val_ind))); % find local minimum in theta_u
ind_sn = val_ind(local_indmin); % map index to full branch

% Convert to a saddle-node point
sn = p_tofold(funcs, branch1.point(ind_sn));

% Get saddle-node calculation method parameters
method = df_mthod(funcs, 'fold');
method.stability.minimal_real_part = -2;

% Correct saddle-node point
[sn, success] = p_correc(funcs, sn, ind_theta_u, [], method.point);
if success ~= 1
    disp('Initial First Saddle-Node Point Correction Unsuccessful')
end

% Compute stability of saddle-node point
sn.stability = p_stabil(funcs, sn, method.stability);

%% Initialise and Continue First Saddle-Node Bifurcation
% Get an empty branch with theta_u and tau as free parameters
branch2 = df_brnch(funcs, [ind_theta_u, ind_taus], 'fold');

% Set bounds for continuation parameters
branch2.parameter.min_bound(1:2,:) = [[ind_theta_u 0.4]' [ind_taus 0]']';
branch2.parameter.max_bound(1:2,:) = [[ind_theta_u 1]' [ind_taus 0.5]']';
branch2.parameter.max_step(1:2,:) = [[ind_theta_u 0.005]' [ind_taus 0.005]']';

% Use saddle-node point as first saddle-node branch point
branch2.point = sn;

% Perturb saddle-node point in tau
sn.parameter(ind_taus) = sn.parameter(ind_taus) + 0.0001;

% Correct saddle-node point in theta_u
[sn, success] = p_correc(funcs, sn, ind_theta_u, [], method.point);
if success ~= 1
    disp('Perturbed First Saddle-Node Point Correction Unsuccessful')
end

% Use perturbed saddle-node point as second branch point
branch2.point(2) = sn;

% Continue hopf branch in both directions
branch2.method.continuation.plot = 0; % hide continuation plot
branch2 = br_contn(funcs, branch2, 500);
branch2 = br_rvers(branch2);
branch2 = br_contn(funcs, branch2, 500);

%% Locate Second Saddle Node Point
% Find the maximum theta_u where us < 0.5
val_ind = find(us < 0.5);  % restrict branch to u < 0.5
[~, local_indmax] = max(arrayfun(@(x) x.parameter(ind_theta_u), ...
    branch1.point(val_ind)));  % find local maximum in theta_u
ind_sn = val_ind(local_indmax); % map index to full branch

% Convert to a saddle-node point
sn = p_tofold(funcs, branch1.point(ind_sn));

% Correct saddle-node point
[sn, success] = p_correc(funcs, sn, ind_theta_u, [], method.point);
if success ~= 1
    disp('Initial Second Saddle-Node Point Correction Unsuccessful')
end

% Compute stability of saddle-node point
sn.stability = p_stabil(funcs, sn, method.stability);

%% Initialise and Continue Second Saddle-Node Bifurcation
% Get an empty branch with theta_u and tau as free parameters
branch3 = df_brnch(funcs, [ind_theta_u, ind_taus], 'fold');

% Set bounds for continuation parameters
branch3.parameter.min_bound(1:2,:) = [[ind_theta_u 0.4]' [ind_taus 0]']';
branch3.parameter.max_bound(1:2,:) = [[ind_theta_u 1]' [ind_taus 0.5]']';
branch3.parameter.max_step(1:2,:) = [[ind_theta_u 0.005]' [ind_taus 0.005]']';

% Use saddle-node point as first saddle-node branch point
branch3.point = sn;

% Perturb saddle-node point in tau
sn.parameter(ind_taus) = sn.parameter(ind_taus) + 0.0001;

% Correct saddle-node point in theta_u
[sn, success] = p_correc(funcs, sn, ind_theta_u, [], method.point);
if success ~= 1
    disp('Perturbed Second Saddle-Node Point Correction Unsuccessful')
end

% Use perturbed saddle-node point as second branch point
branch3.point(2) = sn;

% Continue saddle-node branch in both directions
branch3.method.continuation.plot = 0; % hide continuation plot
branch3 = br_contn(funcs, branch3, 500);
branch3 = br_rvers(branch3);
branch3 = br_contn(funcs, branch3, 500);

%% Locate First Hopf Point
% Find the first hopf bifurcation point
ind_hopf = find(arrayfun(@(x) real(x.stability.l0(1)) > 0, ...
    branch1.point), 1, 'last');

% Convert to a hopf point
hopf = p_tohopf(funcs, branch1.point(ind_hopf));

% Get hopf calculation method parameters
method = df_mthod(funcs, 'hopf');
method.stability.minimal_real_part = -1;

% Correct hopf point
[hopf, success] = p_correc(funcs, hopf, ind_theta_u, [], method.point);
if success ~= 1
    disp('Initial Hopf Point Correction Unsuccessful')
end
first_hopf = hopf; % store hopf point for later use

% Compute stability of hopf point
hopf.stability = p_stabil(funcs, hopf, method.stability);

%% Initialise and Continue Hopf Bifurcation
% Get an empty branch with theta_u and tau as free parameters
branch4 = df_brnch(funcs, [ind_theta_u, ind_taus], 'hopf');

% Set bounds for continuation parameters
branch4.parameter.min_bound(1:2,:) = [[ind_theta_u 0.4]' [ind_taus 0]']';
branch4.parameter.max_bound(1:2,:) = [[ind_theta_u 1]' [ind_taus 0.5]']';
branch4.parameter.max_step(1:2,:) = [[ind_theta_u 0.005]' [ind_taus 0.005]']';

% Use hopf point as first hopf branch point
branch4.point = hopf;

% Perturb hopf point in theta_u
hopf.parameter(ind_theta_u) = hopf.parameter(ind_theta_u) - 0.0001;

% Correct hopf point in tau
[hopf, success] = p_correc(funcs, hopf, ind_taus, [], method.point);
if success ~= 1
    disp('Perturbed Hopf Point Correction Unsuccessful')
end

% Use perturbed hopf point as second branch point
branch4.point(2) = hopf;

% Continue hopf branch in both directions
branch4.method.continuation.plot = 0; % hide continuation plot
branch4 = br_contn(funcs, branch4, 500);
branch4 = br_rvers(branch4);
branch4 = br_contn(funcs, branch4, 500);

%% Initialise and Continue Periodic Orbits
% Construct small amplitude oscillation around first hopf point
intervals = 18;
degree = 3;
[psol, stepcond] = p_topsol(funcs, first_hopf, 1e-2, degree, intervals);

% Correct periodic solution
method = df_mthod(funcs, 'psol');
[psol, success] = p_correc(funcs, psol, ind_theta_u, stepcond, method.point);
if success ~= 1
    disp('Periodic Solution Correction Unsuccessful')
end

% Get an empty branch with theta_u as a free parameter
branch5 = df_brnch(funcs, ind_theta_u, 'psol');

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

% Continue periodic solutions branch
branch5.method.continuation.plot = 0; % hide continuation plot
branch5 = br_contn(funcs, branch5, 500);

%% Initialise Folds of Periodic Orbits
% Speed up computations by vectorisation
sys_rhs = @(xx, par) [ ...
    -xx(1,1,:) + 1 / (1+ exp(-par(2) * (par(7) + par(3) * xx(1,2,:) + par(4) * xx(2,2,:))));....
    par(1) * (-xx(2,1,:) + 1 / (1+ exp(-par(2) * (par(8) + par(5) * xx(1,2,:) + par(6) * xx(2,2,:)))))];
vfuncs = set_funcs( ...
    'sys_rhs', sys_rhs, ...
    'sys_tau', @() 9, ...
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

%% Initialise and Continue First Saddle-Node of Periodic Orbits
psol = branch5.point(indmax:indmax+1);
intervals = 40; degree = 4;
psol = arrayfun(@(p) p_remesh(p, degree, intervals), psol); % refine
method.point.adapt_mesh_after_correct = 1;
method.point.newton_max_iterations = 7;
method.point.newton_nmon_iterations = 2;
psol = arrayfun(@(p) p_correc(funcs, p, [], [], method.point), psol); % correct

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

% Continue first saddle-node
branch6.method.continuation.plot = 0; % hide continuation plot
branch6.method.point.print_residual_info = 0;
branch6 = br_contn(foldfuncs, branch6, 100);
branch6 = br_rvers(branch6);
branch6 = br_contn(foldfuncs, branch6, 100);

%% Initialise and Continue Second Saddle-Node of Periodic Orbits
psol = branch5.point(indmin-1:indmin);
intervals = 40; degree = 4;
psol = arrayfun(@(p) p_remesh(p, degree, intervals), psol); % refine
method.point.adapt_mesh_after_correct = 1;
method.point.newton_max_iterations = 7;
method.point.newton_nmon_iterations = 2;
psol = arrayfun(@(p) p_correc(funcs, p, [], [], method.point), psol); % correct

% Copy branch 5 and start at min theta_u point
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

% Continue second saddle-node
branch7.method.continuation.plot = 0; % hide continuation plot
branch7.method.point.print_residual_info = 0;
branch7 = br_contn(foldfuncs, branch7, 100);
branch7 = br_rvers(branch7);
branch7 = br_contn(foldfuncs, branch7, 100);

%% Construct Final Figure
% Get deafult plotting measures: x = theta_u, y = tau
[xm, ym] = df_measr(0, branch2);

% Extract data from branches
    % First saddle-node bifurcation
    bifn.sn1.x = br_measr(branch2, xm);
    bifn.sn1.y = br_measr(branch2, ym);

    % Second saddle-node bifurcation
    bifn.sn2.x = br_measr(branch3, xm);
    bifn.sn2.y = br_measr(branch3, ym);

    % Hopf bifurcation
    bifn.hopf.x = br_measr(branch4, xm);
    bifn.hopf.y = br_measr(branch4, ym);

    % First saddle-node of periodic orbits
    bifn.snpo1.x = br_measr(branch6, xm);
        % Regularise and restrict to crossing point
        bifn.snpo1.x = linspace(max(bifn.snpo1.x), min(bifn.snpo1.x), 19);
        bifn.snpo1.x = bifn.snpo1.x(1:12);
    bifn.snpo1.y = br_measr(branch6, ym);
        % Regularise and restrict to crossing point
        bifn.snpo1.y = linspace(min(bifn.snpo1.y), max(bifn.snpo1.y), 19);
        bifn.snpo1.y = bifn.snpo1.y(1:12);

    % Second saddle-node of periodic orbits
    bifn.snpo2.x = br_measr(branch7, xm);
        % Regularise and restrict to crossing point
        bifn.snpo2.x = linspace(min(bifn.snpo2.x), max(bifn.snpo2.x), 19);
        bifn.snpo2.x = bifn.snpo2.x(1:12);
    bifn.snpo2.y = br_measr(branch7, ym);
        % Regularise and restrict to crossing point
        bifn.snpo2.y = linspace(min(bifn.snpo2.y), max(bifn.snpo2.y), 19);
        bifn.snpo2.y = bifn.snpo2.y(1:12);
end