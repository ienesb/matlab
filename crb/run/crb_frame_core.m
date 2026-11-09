function res = crb_frame_core(ch, gcp, grcp, o)
% CRB_FRAME_CORE  Monte-Carlo genie-CRB sweep shared by run_crb_eva and
% run_crb_mt40_frame: for nT trials draw a channel realization (chan_realize)
% and a symbol grid, build the FIM for every waveform / ICI case at sigma2 = 1,
% and store per-path CRBs, condition numbers, validity, the closed-form
% (E J)^-1 for the CP waveforms and the continuous-time floor.
%   o fields: nT, waveforms (cell), ici_list, modtype, seed, kappa_thresh
% CP waveforms use gcp, RCP waveforms grcp (same N, M, deltaf; TS differs).
    N = gcp.N;  M = gcp.M;  P = ch.P;
    nW = numel(o.waveforms);  nI = numel(o.ici_list);  nT = o.nT;
    crb_tau = nan(nT, P, nW, nI);  crb_nu = crb_tau;  kappa = nan(nT, nW, nI);  valid = false(nT, nW, nI);
    exp_tau = nan(nT, P, nW, nI);  exp_nu = exp_tau;
    fl_tau = nan(nT, P);  fl_nu = fl_tau;  alphas = nan(nT, P);  nus = nan(nT, P);
    wfs = o.waveforms;  ici_list = o.ici_list;  modtype = o.modtype;  seed = o.seed;  kth = o.kappa_thresh;
    t0 = tic;
    parfor t = 1:nT
        rs = crn_stream(seed, 0, t);
        cr = chan_realize(ch, gcp, rs);
        S  = gen_symbols(N, M, modtype, rs);
        ct = ct_crlb_floor(gcp, cr, 1);
        fl_tau(t, :) = ct.tau_floor.';  fl_nu(t, :) = ct.nu_floor.';
        alphas(t, :) = cr.alpha.';  nus(t, :) = cr.nu.';
        ctau = nan(P, nW, nI);  cnu = ctau;  kap = nan(nW, nI);  val = false(nW, nI);  etau = ctau;  enu = ctau;
        for iw = 1:nW
            wf = wfs{iw};
            if startsWith(wf, 'cp'), g = gcp; else, g = grcp; end
            for ii = 1:nI
                ici = ici_list(ii);
                J = fim(S, wf, g, cr, ici, 1);
                Sc = crb_from_fim(J, g, cr);
                ctau(:, iw, ii) = Sc.crb_tau;  cnu(:, iw, ii) = Sc.crb_nu;
                kap(iw, ii) = Sc.kappa;  val(iw, ii) = Sc.valid && Sc.kappa <= kth;
                if startsWith(wf, 'cp')
                    Se = crb_from_fim(fim_expected(wf, g, cr, ici, 1), g, cr);
                    etau(:, iw, ii) = Se.crb_tau;  enu(:, iw, ii) = Se.crb_nu;
                end
            end
        end
        crb_tau(t, :, :, :) = ctau;  crb_nu(t, :, :, :) = cnu;  kappa(t, :, :) = kap;  valid(t, :, :) = val;
        exp_tau(t, :, :, :) = etau;  exp_nu(t, :, :, :) = enu;
    end
    res = struct('crb_tau', crb_tau, 'crb_nu', crb_nu, 'exp_tau', exp_tau, 'exp_nu', exp_nu, ...
        'kappa', kappa, 'valid', valid, 'fl_tau', fl_tau, 'fl_nu', fl_nu, 'alphas', alphas, 'nus', nus, ...
        'waveforms', {o.waveforms}, 'ici_list', o.ici_list, 'nT', nT, 'gcp', gcp, 'grcp', grcp, 'ch', ch, ...
        'seconds', toc(t0));
    res.drop_frac = reshape(1 - mean(valid, 1), nW, nI);
    fprintf('  core: %d trials in %.0f s (%.1f s/trial)\n', nT, res.seconds, res.seconds/nT);
end
