function cr = chan_realize(ch, g, rs)
% CHAN_REALIZE  Draw one channel realization from a chan_* specification.
%   cr = chan_realize(ch, g, rs)      rs: RandStream (crn_stream) or [] = global
% Output cr: alpha (P x 1 complex), tau [s], nu [Hz], P, eta_bar (4P x 1),
%            l_tap, k_real = nu*M*TS, name.
% EVA:  alpha ~ CN(0,P_l) ('rayleigh') or sqrt(P_l) e^{j phi} ('constmag');
%       nu_l = doppler_scale * nu_max * cos(theta_l), theta_l ~ U[-pi,pi);
%       legacy preset: deterministic (fixed theta, phi).
% MT40: fixed geometry and nu; gains 1 / e^{j phi} / CN(0,1).
% Controls: one tap at tau = nu = 0.
    if nargin < 3, rs = []; end
    P = ch.P;
    switch ch.family
        case 'eva'
            if ch.legacy
                alpha = 10.^(ch.pow_db/20) .* exp(1j*ch.phi);
                nu    = ch.nu_max * cos(ch.theta);
            else
                alpha = draw_gains(ch.gain_model, ch.P_lin, rs);
                theta = 2*pi*urand(P, rs) - pi;
                nu    = ch.doppler_scale * ch.nu_max * cos(theta);
            end
            tau = ch.tau;
            if ch.merge_colliding && isfield(ch, 'merge_groups')
                for i = 1:numel(ch.merge_groups)
                    grp = ch.merge_groups{i};
                    alpha(grp(1)) = sum(alpha(grp));   % coherent sum on the shared tap
                    nu(grp(1))    = mean(nu(grp));     % documented approximation
                end
                alpha = alpha(ch.merge_keep);  tau = tau(ch.merge_keep);  nu = nu(ch.merge_keep);
                P = numel(alpha);
            end
        case 'mt40'
            switch ch.gain_model
                case 'unit',      alpha = ones(P, 1);
                case 'randphase', alpha = exp(1j*2*pi*urand(P, rs));
                case 'cn01',      alpha = draw_gains('rayleigh', ones(P,1), rs);
                otherwise, error('chan_realize:gain', 'unknown gain model %s', ch.gain_model);
            end
            tau = ch.tau;  nu = ch.nu;
        case 'control'
            switch ch.kind
                case 'awgn',          alpha = exp(1j*2*pi*urand(1, rs));
                case 'flat_rayleigh', alpha = draw_gains('rayleigh', 1, rs);
            end
            tau = 0;  nu = 0;
        otherwise
            error('chan_realize:family', 'unknown channel family %s', ch.family);
    end
    cr.name  = ch.name;
    cr.P     = P;
    cr.alpha = alpha(:);
    cr.tau   = tau(:);
    cr.nu    = nu(:);
    cr.eta_bar = eta_pack(cr.alpha, cr.tau, cr.nu);
    cr.l_tap  = round(cr.tau / g.Ts);
    cr.k_real = cr.nu * g.M * g.TS;
end

function a = draw_gains(model, P_lin, rs)
    P_lin = P_lin(:);  P = numel(P_lin);
    switch lower(model)
        case 'rayleigh'
            a = sqrt(P_lin/2) .* (nrand(P, rs) + 1j*nrand(P, rs));
        case 'constmag'
            a = sqrt(P_lin) .* exp(1j*2*pi*urand(P, rs));
        otherwise
            error('chan_realize:gain', 'unknown gain model %s', model);
    end
end

function u = urand(n, rs)
    if isempty(rs), u = rand(n, 1); else, u = rand(rs, n, 1); end
end

function z = nrand(n, rs)
    if isempty(rs), z = randn(n, 1); else, z = randn(rs, n, 1); end
end
