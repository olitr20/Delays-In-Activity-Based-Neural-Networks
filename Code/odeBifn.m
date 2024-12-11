function [hopf, sn, bt] = odeBifn(hopf, sn, bt, p)
% BIFN_ODE  Calculate hopf, saddle-node, and bogdanov-takens bifurcations
% in the \theta_{u}, \theta_{v} plane.
%   Input:
%       hopf: structure containing us - range of u values to search for
%           hopf bifurcations. Found by analysing Tr L (eq. 2.5) in the u-v
%           plane, taking absolute lower and upper bounds as limits of Tr L
%           in the real plane and relative upper and lower bounds as roots
%           of the Tr L.
%       sn: structure containing us - range of u values to search for
%           saddle-node bifurcations. Found by analysing Det L, calculated
%           from L (eq. 2.4), in the u-v plane, taking absolute lower and
%           upper bounds as the roots of Det L. Note that Det L = \alpha -
%           \alpha d \beta v^{*} (1 - v^{*}) - \alpha a \beta u^{*}
%           (1 - u^{*}) + \alpha \beta^{2} (ad - bc) u^{*} (1 - u^{*})
%           v^{*} (1 - v^{*}).
%       bt: structure containing us - range of u values to searcg for
%           bogdanov-takens bifurcations. Found by analysing Tr L = 0 and
%           Det L = 0
%       p:  structure containing model parameters of the form {\alpha,
%           \beta, a, b, c, d}
%   Output:
%       hopf: structure containing us, v_plus, v_minus, theta. Theta is a
%           structure containing uP (u plus), uM (u minus), vP (v plus) and
%           vM (v minus).
%       sn: structure containing us, v_plus, v_minus, theta. Theta is a
%           structure containing uP (u plus), uM (u minus), vP (v plus) and
%           vM (v minus).
%       bt: structure containing us, v_plus, v_minus, theta. Theta is a
%           structure containing uP (u plus), uM (u minus), vP (v plus) and
%           vM (v minus).

    switch nargin
        case 3 % if only hopf, sn and p are provided
            p = bt;
            bt = [];
        case 4
            % all inputs provided
    end

    % Calculate hopf bifurcation
    [hopf.v_plus, hopf.v_minus] = vTrace(hopf.us,p);
    hopf.theta = fixedPoint(hopf.us,p,hopf.v_plus,hopf.v_minus);
    hopf.theta.uP = real(hopf.theta.uP); hopf.theta.uM = real(hopf.theta.uM);
    hopf.theta.vP = real(hopf.theta.vP); hopf.theta.vM = real(hopf.theta.vM);
    
    % Calculate saddle node bifurcation
    [sn.v_plus, sn.v_minus] = vDet(sn.us,p);
    sn.theta = fixedPoint(sn.us,p,sn.v_plus,sn.v_minus);
    sn.theta.uP = real(sn.theta.uP); sn.theta.uM = real(sn.theta.uM);
    sn.theta.vP = real(sn.theta.vP); sn.theta.vM = real(sn.theta.vM);
    
    
    % Calculate bogdanov-takens bifurcations
    if ~isempty(bt)
        [bt.v_plus, bt.v_minus] = vTrace(bt.us,p);
        bt.theta = fixedPoint(bt.us, p, bt.v_plus,bt.v_minus);
        bt.theta.uP = real(bt.theta.uP); bt.theta.uM = real(bt.theta.uM);
        bt.theta.vP = real(bt.theta.vP); bt.theta.vM = real(bt.theta.vM);
    end

%% --------------------------------------------------------------------- %%
% ----------------------------- f_inv(x,p) ------------------------------ %
    % Define the inverse sigmoid function, for z = u,v
    function f = f_inv(z,p)
        f = (1 ./ p.beta) .* log(z ./ (1 - z));
    end

% ---------------------- Hopf Fixed Point Equation ---------------------- %
    % Define the fixed point equation derived from Tr L
    function [v_trace_plus, v_trace_minus] = vTrace(us,p)
        sqrt_term = 1 - (4 .* ...
            (((1 + p.alpha) ./ p.beta) - p.a .* us .* (1 - us))) ...
            ./ (p.alpha .* p.d);
        v_trace_plus = (1 + sqrt(sqrt_term)) ./ 2;
        v_trace_minus = (1 - sqrt(sqrt_term)) ./ 2;
    end

% ------------------- Saddle Node Fixed Point Equation ------------------ %
    % Define the fixed point equation derived from Det L
    function [v_det_plus, v_det_minus] = vDet(us,p)
        gamma = (p.d .* p.beta) - (p.beta .^ 2) .* ...
            (p.a .* p.d - p.b * p.c) .* us .* (1 - us);
        delta = 1 - (p.a .* p.beta .* us .* (1 - us));
        sqrt_term = (gamma .^ 2) - 4 .* gamma .* delta;

        v_det_plus = (gamma + sqrt(sqrt_term)) ./ (2 .* gamma);
        v_det_minus = (gamma - sqrt(sqrt_term)) ./ (2 .* gamma);
    end

% ----------------------- Parameterised Equations ----------------------- %
    % Define the pair of equations for θ_u and θ_v
    function theta = fixedPoint(us,p,vPlus,vMinus)
        theta.uP = f_inv(us,p) - (p.a .* us) - (p.b .* vPlus);
        theta.uM = f_inv(us,p) - (p.a .* us) - (p.b .* vMinus);
        theta.vP = f_inv(vPlus,p) - (p.c .* us) - (p.d .* vPlus);
        theta.vM = f_inv(vMinus,p) - (p.c .* us) - (p.d .* vMinus);
    end
end