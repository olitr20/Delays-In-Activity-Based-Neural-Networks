%% Initialize continuation of Hopf bifurcations
%%
function [hbranch,suc]=SetupHopf(funcs,branch,ind,varargin)
%% Inputs
%
% * |funcs|: functions used for DDE
% * |branch|: branch of stst along which Hopf was discovered
% * |ind|: index of approximate Hopf point
%
% Important name-value pair inputs
%
% * |'contpar'| (integers default |[]|): index of continuation parameters
%  (replacing free pars of branch)
% * |'correc'| (logical, default true): apply |p_correc| to first points on
% hopf branch
% * |'dir'| (integer, default |[]|): which parameter to vary initially
% along Hopf branch (hbranch has only single point if dir is empty)
% * |'step'| (real, default |1e-3|): size of initial step if dir is non-empty
%
% All other named arguments are passed on to fields of |hbranch|
%% Outputs
% 
% * |hbranch|: Hopf branch with first point (or two points)
% * |suc|: flag whether corection was successful
%
% Parameter limits for etc are inherited from branch, unless overridden by
% optional input arguments.
%
% (c) DDE-BIFTOOL v. 3.1.1(56), 29/06/2014
%
%% wrapper around SetupStstBifurcation to ensure backward compatibility
[hbranch,suc]=SetupStstBifurcation(funcs,branch,ind,'hopf',varargin{:});
end
