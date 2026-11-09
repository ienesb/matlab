function S = crb_from_fim(J, g, cr, varargin)
% CRB_FROM_FIM  Conditioned CRB extraction from a 4P x 4P FIM.
%   S = crb_from_fim(J, g, cr, 'name', value, ...)
% Options (defaults):
%   'active'   auto     logical 4P x 1 mask of estimated parameters; auto = all,
%                       minus the nu entries when every cr.nu == 0 (structurally
%                       null Doppler block, e.g. MT40)
%   'keep_nu'  false    force the nu block active even if nu == 0
%   'method'   'bin'    'bin': exact change to bin units J_bin = D J D
%                       (eta_scale) then pinv, back to SI (default);
%                       'jacobi': legacy pinv_equilibrated on the SI matrix
%   'efi'      false    Schur-complement the gain block first (equivalent
%                       Fisher information); same numbers as the block of the
%                       full inverse when J is full rank
% Output S: crb_tau, crb_nu (P x 1, NaN when inactive), C (4P x 4P CRB in SI,
% zero rows/cols for inactive), kappa = cond of the active bin-unit FIM,
% active, valid (finite and PSD), tr_tau, tr_nu.
    o = struct('active', [], 'keep_nu', false, 'method', 'bin', 'efi', false);
    for i = 1:2:numel(varargin), o.(varargin{i}) = varargin{i+1}; end
    P = numel(cr.alpha);
    idx = eta_index(P);
    n = 4*P;
    if isempty(o.active)
        active = true(n, 1);
        if all(cr.nu == 0) && ~o.keep_nu, active(idx.nu) = false; end
    else
        active = logical(o.active(:));
    end
    J = (J + J.')/2;
    C = zeros(n);
    switch lower(o.method)
        case 'bin'
            D  = eta_scale(g, P);
            Jb = D*J*D;
            Ja = Jb(active, active);
            if o.efi
                ia = ismember(find(active), idx.gain);          % gains among active
                Jaa = Ja(ia, ia);  Jak = Ja(ia, ~ia);  Jkk = Ja(~ia, ~ia);
                Jschur = Jkk - Jak' * (Jaa \ Jak);
                Ck = pinv((Jschur + Jschur')/2);
                Ca = zeros(size(Ja));  Ca(~ia, ~ia) = Ck;        % gain CRBs not returned in EFI mode
            else
                Ca = pinv((Ja + Ja')/2);
            end
            Cb = zeros(n);  Cb(active, active) = Ca;
            C  = D*Cb*D;
            kappa = cond(Ja);
        case 'jacobi'
            Ja = J(active, active);
            d = sqrt(abs(diag(Ja)));  d(d <= 0) = 1;
            Dinv = diag(1./d);
            Ca = Dinv * pinv(Dinv*Ja*Dinv) * Dinv;
            C(active, active) = Ca;
            kappa = cond(Dinv*Ja*Dinv);
        otherwise
            error('crb_from_fim:method', 'unknown method %s', o.method);
    end
    dC = diag(C);
    S.crb_tau = dC(idx.tau);  S.crb_tau(~active(idx.tau)) = NaN;
    S.crb_nu  = dC(idx.nu);   S.crb_nu(~active(idx.nu))   = NaN;
    S.C = C;  S.kappa = kappa;  S.active = active;
    S.tr_tau = sum(S.crb_tau, 'omitnan');
    S.tr_nu  = sum(S.crb_nu,  'omitnan');
    S.valid  = all(isfinite(dC(active))) && all(dC(active) >= -1e-12*max(abs(dC(active))));
end
