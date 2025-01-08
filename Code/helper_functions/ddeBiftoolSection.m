function [stst, po] = ddeBiftoolSection(p)
% DDEBIFTOOLSECTION  Calculate steady state and periodic orbit branches
% in a two parameter space.
%   Input:
%       p:  structure containing model parameters of the form {\alpha,
%           \beta, a, b, c, d, \theta_{u}, \theta_{v}, \tau_{1}, \tau_{2}}.
%   Output:
%       stst: structure containing kind (type of point, 'stst'), parameter
%           (parameter values used for continuation), x (initial guess for
%           steady state location), stability (calculated stability of
%           every point along steady state branch), stable1 (structure
%           containing x and y values of the first stable section of the
%           steady state branch in the two parameter space), unstable
%           (structure containing x and y values of the unstable section of
%           the steady state branch in the two parameter space), stable2
%           (structure containing x and y values of the second stable
%           section of the steady state branch in the two parameter space.
%       po: structure containing stable1 (structure containing x and y
%           values of the first stable section of the periodic orbit branch
%           in the two parameter space), unstable (structure containing x
%           and y values of the unstable section of the periodic orbit
%           branch in the two parameter space), stable2 (structure
%           containing x and y values of the second stable section of the
%           periodic orbit branch in the two parameter space.

%% Initislise System
addpath('helper_functions/ddebiftool/');
addpath('helper_functions/ddebiftool_utilities/');
addpath('helper_functions/ddebiftool_extra_psol/');

%#ok<*ASGLU,*NOPTS,*NASGU>

%% Define Model
% Right-hand side
sys_rhs = @(xx, par) [ ...
    -xx(1,1) + 1 / (1 + exp(-par(2) * (par(7) + par(3) * xx(1,2) + par(4) * xx(2,2)))); ...
    par(1) * (-xx(2,1) + 1 / (1 + exp(-par(2) * (par(8) + par(5) * xx(1,2) + par(6) * xx(2,2)))))];

% Delays and Continuation Parameters
sys_tau = @() 9;
ind_theta_u = 7;
ind_taus = 9;

% Define system functions
funcs = set_funcs( ...
    'sys_rhs', sys_rhs, ...
    'sys_tau', sys_tau);

% Define continuation parameters
if p.tau_1 == 0.09; min_theta_u = 0.61; max_theta_u = 1;
elseif p.tau_1 == 0.2; min_theta_u = 0.4; max_theta_u = 0.85;
elseif p.tau_1 == 0.5; min_theta_u = 0.4; max_theta_u = 0.97;
else; min_theta_u = 0.4; max_theta_u = 1;
end

%% Initial Guess for Steady State
% Initialise steady state point
stst.kind='stst';
stst.parameter=[p.alpha, p.beta, p.a, p.b, p.c, p.d, ...
    p.theta_u, p.theta_v, p.tau_1];

% Extract u and v values from odeSim endpoint
[stst_u, stst_v] = odeSim(p);
stst.x=[stst_u; stst_v];

method = df_mthod(funcs, 'stst');
method.stability.minimal_real_part = -2;

% Correct steady state point
[stst, success] = p_correc(funcs, stst, [], [], method.point);
if success ~= 1
    disp('Initial Steady State Correction Unsuccessful')
end

% Compute point stability
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
if success ~= 1
    disp('Perturbed Steady State Correction Unsuccessful')
end

% Use new steady state point as second steady state branch point
branch1.point(2) = stst;

% Continue steady state branch in both directions
branch1.method.continuation.plot = 0; % hide continuation plot
branch1 = br_contn(funcs, branch1, 200);
branch1 = br_rvers(branch1);
branch1 = br_contn(funcs, branch1, 200);

% Calculate stability along branch
branch1.method.stability.minimal_real_part = -2;
branch1 = br_stabl(funcs, branch1, 0, 0);

%% Locate the First Hopf Point
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
if success ~= 1
    disp('Initial Hopf Point Correction unsuccessful')
end
first_hopf = hopf; % store hopf point for later use

%% Initialise and Continue Periodic Orbits
% Construct small amplitude oscillation
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
branch2 = df_brnch(funcs, ind_theta_u, 'psol');

% Set bounds for continuation parameter
branch2.parameter.min_bound(1,:) = [ind_theta_u min_theta_u];
branch2.parameter.max_bound(1,:) = [ind_theta_u max_theta_u];
branch2.parameter.max_step(1,:) = [ind_theta_u 0.005];

% Make degenerate periodic solution with amplitude zero at hopf point
deg_psol = p_topsol(funcs, first_hopf, 0, degree, intervals);

% Use deg_psol and psol as first two points on branch
deg_psol.mesh = []; % clear the mesh field to save memory and avoid adaptive mesh selection
branch2.point = deg_psol;
psol.mesh = [];
branch2.point(2) = psol;

% Continue periodic solutions branch
branch2.method.continuation.plot = 0; % hide continuation plot
branch2 = br_contn(funcs, branch2, 200);

%% Extract max u and theta_u along the branch
% Get default measures for branch 2
[xm, ym] = df_measr(0, branch2);

% Update measures
ym.field = 'profile'; ym.row = 1; ym.col = 'all';

% Extract periodic orbits branch
po_x = br_measr(branch2,xm); % theta_us
po_y = br_measr(branch2,ym); % us
po_y = max(po_y, [], 2); % find maximum u

% Calculate stability along branch
branch2.method.stability.minimal_real_part = -2;
branch2 = br_stabl(funcs, branch2, 0, 0);

%% Extract Data for Final Figure
% Get default measures for branch 1
[xm, ym] = df_measr(0, branch1);

% Extract steady states branch
stst_x = br_measr(branch1, xm);
stst_y = br_measr(branch1, ym);

% Locate first hopf bifurcation
ind_hopf_1 = find(arrayfun(@(x) real(x.stability.l0(1)) > 0, ...
    branch1.point), 1, 'first');

% Locate second hopf bifurcation
ind_hopf_2 = find(arrayfun(@(x) real(x.stability.l0(1)) > 0, ...
    branch1.point), 1, 'last');

% Extract first stable section of steady state branch
stst.stable1.x = stst_x(1:ind_hopf_1);
stst.stable1.y = stst_y(1:ind_hopf_1);

% Extract unstable section of steady state branch
stst.unstable.x = stst_x(ind_hopf_1:ind_hopf_2);
stst.unstable.y = stst_y(ind_hopf_1:ind_hopf_2);

% Extract second stable section of steady state branch
stst.stable2.x = stst_x(ind_hopf_2:end);
stst.stable2.y = stst_y(ind_hopf_2:end);

% Get default measures for branch 2
[xm, ym] = df_measr(0, branch2);

% Update measures
ym.field = 'profile'; ym.row = 1; ym.col = 'all';

% Extract periodic orbits branch
po_x = br_measr(branch2,xm); % theta_us
po_y = br_measr(branch2,ym); % us
po_y = max(po_y, [], 2); % find maximum u

% Locate first hopf bifurcation of periodic orbits
ind_hopf_po_1 = find(arrayfun(@(x) real(x.stability.mu(1)) > 1.1, ...
    branch2.point), 1, 'first');

% Locate second hopf bifurcation of periodic orbits
ind_hopf_po_2 = find(arrayfun(@(x) real(x.stability.mu(1)) > 1.1, ...
    branch2.point), 1, 'last');

% Account for no instability (no hopf bifurcations of periodic orbits)
if isempty(ind_hopf_po_1)
    ind_hopf_po_1 = length(branch2.point);
end

% Extract first stable section of periodic orbit branch
if ind_hopf_po_1 ~= 1 % account for no stability
    po_x_stable_1 = po_x(1:ind_hopf_po_1);
    po_y_stable_1 = po_y(1:ind_hopf_po_1);
    [po.stable1.x, po.stable1.y] = regularise(po_x_stable_1, po_y_stable_1);
end

% Extract unstable section of periodic orbit branch
if ind_hopf_po_1 ~= length(branch2.point) % account for no instability
    po_x_unstable = po_x(ind_hopf_po_1:ind_hopf_po_2);
    po_y_unstable = po_y(ind_hopf_po_1:ind_hopf_po_2);
    [po.unstable.x, po.unstable.y] = regularise(po_x_unstable, po_y_unstable);
end

% Extract second stable section of periodic orbit branch
if ind_hopf_po_1 ~= length(branch2.point) % account for no instability
    if ind_hopf_po_1 ~= 1 % account for a single stable branch already extracted
        po_x_stable_2 = po_x(ind_hopf_po_2:end);
        po_y_stable_2 = po_y(ind_hopf_po_2:end);
        [po.stable2.x, po.stable2.y] = regularise(po_x_stable_2, po_y_stable_2);
    end
end

%% --------------------------------------------------------------------- %%
% --------------------------- regularise(x,y) --------------------------- %
    % Define the inverse sigmoid function, for z = u,v
    function [reg_x, reg_y] = regularise(x,y)
        % Calculate cumulative euclidean distances
        dist = sqrt(diff(x).^2 + diff(y).^2);
        cum_dist = [0; cumsum(dist)];

        % Define uniform spacing
        tot_dist = cum_dist(end);
        N = round(length(x) * 0.3);
        reg_dist = linspace(0, tot_dist, N);

        % Interpolate onto uniform spacing
        reg_x = interp1(cum_dist, x, reg_dist, 'linear');
        reg_y = interp1(cum_dist, y, reg_dist, 'linear');
    end
end