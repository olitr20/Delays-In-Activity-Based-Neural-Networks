function [hopf, sn, bt] = odeBifn(p)
% ODEBIFN  Calculate hopf, saddle-node, and bogdanov-takens bifurcations
% in the \theta_{u}, \theta_{v} plane.
%   Input:
%       p:  structure containing model parameters of the form {\alpha,
%           \beta, a, b, c, d}.
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

    % Define u ranges for bifurcation search
    hopf.us = 0:0.00001:1;
    sn.us = 0:0.00001:1;

    % Extract only valid u values
    [hopf.v_plus, hopf.v_minus] = vTrace(hopf.us,p);
    hopf.us = hopf.us( ...
        imag(hopf.v_plus) == 0 & real(hopf.v_plus) < 1 & ...
        imag(hopf.v_minus) == 0 & real(hopf.v_minus) > 0);

    [sn.v_plus, sn.v_minus] = vDet(sn.us,p);
    sn.us = sn.us( ...
        imag(sn.v_plus) == 0 & real(sn.v_plus) > 0);

    % Calculate hopf bifurcation
    [hopf.v_plus, hopf.v_minus] = vTrace(hopf.us,p);
    hopf.theta = fixedPoint(hopf.us,p,hopf.v_plus,hopf.v_minus);
    hopf.theta.uP = real(hopf.theta.uP); hopf.theta.uM = real(hopf.theta.uM);
    hopf.theta.vP = real(hopf.theta.vP); hopf.theta.vM = real(hopf.theta.vM);

    % Separete non-contiguous sections
    jumps = abs(diff(hopf.theta.uP)) > 100 * mean(abs(diff(hopf.theta.uP)));
    hopf.theta.uP([false, jumps]) = NaN;

    jumps = abs(diff(hopf.theta.uM)) > 100 * mean(abs(diff(hopf.theta.uM)));
    hopf.theta.uM([false, jumps]) = NaN;
    
    % Calculate saddle node bifurcation
    [sn.v_plus, sn.v_minus] = vDet(sn.us,p);
    sn.theta = fixedPoint(sn.us,p,sn.v_plus,sn.v_minus);
    sn.theta.uP = real(sn.theta.uP); sn.theta.uM = real(sn.theta.uM);
    sn.theta.vP = real(sn.theta.vP); sn.theta.vM = real(sn.theta.vM);

    % Separete non-contiguous sections
    jumps = abs(diff(sn.theta.uP)) > 100 * mean(abs(diff(sn.theta.uP)));
    sn.theta.uP([false, jumps]) = NaN;

    jumps = abs(diff(sn.theta.uM)) > 100 * mean(abs(diff(sn.theta.uM)));
    sn.theta.uM([false, jumps]) = NaN;
    
    % Locate bogdanov-takens bifurcations
    try
        options = optimoptions('fsolve','Display','none','Algorithm','levenberg-marquardt');
        bt.us(1) = fsolve(@uTraceDet, 0, options);
        bt.us(2) = fsolve(@uTraceDet, 1, options);

        % Calculate bogdanov-takens bifurcations
        [bt.v_plus, bt.v_minus] = vTrace(bt.us,p);
        bt.theta = fixedPoint(bt.us, p, bt.v_plus,bt.v_minus);
        bt.theta.uP = real(bt.theta.uP); bt.theta.uM = real(bt.theta.uM);
        bt.theta.vP = real(bt.theta.vP); bt.theta.vM = real(bt.theta.vM);
    catch
        bt = [];
    end

%% --------------------------------------------------------------------- %%
% ----------------------------- f_inv(z,p) ------------------------------ %
    % Define the inverse sigmoid function, for z = u,v
    function f = f_inv(z,p)
        f = (1 ./ p.beta) .* log(z ./ (1 - z));
    end

% ---------------------------- vTrace(us,p) ----------------------------- %
    % Define the conditions for a hopf bifurcation
    function [v_trace_plus, v_trace_minus] = vTrace(us,p)
        sqrt_term = 1 - (4 .* ...
            (((1 + p.alpha) ./ p.beta) - p.a .* us .* (1 - us))) ...
            ./ (p.alpha .* p.d);
        v_trace_plus = (1 + sqrt(sqrt_term)) ./ 2;
        v_trace_minus = (1 - sqrt(sqrt_term)) ./ 2;
    end

% ----------------------------- vDet(us,p) ------------------------------ %
    % Define the conditions for a saddle-node bifurcation
    function [v_det_plus, v_det_minus] = vDet(us,p)
        gamma = (p.d .* p.beta) - (p.beta .^ 2) .* ...
            (p.a .* p.d - p.b * p.c) .* us .* (1 - us);
        delta = 1 - (p.a .* p.beta .* us .* (1 - us));
        sqrt_term = (gamma .^ 2) - 4 .* gamma .* delta;

        v_det_plus = (gamma + sqrt(sqrt_term)) ./ (2 .* gamma);
        v_det_minus = (gamma - sqrt(sqrt_term)) ./ (2 .* gamma);
    end

% ---------------------------- uTraceDet(us) ---------------------------- %
    % Define the conditions for a bogdanov-takens bifurcation
    function bt_us = uTraceDet(us)
        [v_trace_plus, v_trace_minus] = vTrace(us, p);
        [v_det_plus, v_det_minus] = vDet(us, p);

        bt_us(1) = v_trace_plus - v_det_minus;  % Match one branch
        bt_us(2) = v_trace_minus - v_det_plus;  % Match the other branch
    end

% ------------------- fixed-Point(us,p,vPlus,vMinus) -------------------- %
    % Define the fixed point equations for \theta_{u} and \theta_{v}
    function theta = fixedPoint(us,p,vPlus,vMinus)
        theta.uP = f_inv(us,p) - (p.a .* us) - (p.b .* vPlus);
        theta.uM = f_inv(us,p) - (p.a .* us) - (p.b .* vMinus);
        theta.vP = f_inv(vPlus,p) - (p.c .* us) - (p.d .* vPlus);
        theta.vM = f_inv(vMinus,p) - (p.c .* us) - (p.d .* vMinus);
    end
end