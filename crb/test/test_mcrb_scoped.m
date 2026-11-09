function r = test_mcrb_scoped()
% Independent identities for publication validation, without old/impl oracles.
    r = [];
    g = make_grid(32,16,15e3,4/(32*15e3),'cp');
    cr = struct('P',2,'alpha',[1;.6*exp(.7j)], ...
        'tau',[.4;3.2]*g.Ts,'nu',[.2;-.35]*g.bin_nu);
    cr.eta_bar = eta_pack(cr.alpha,cr.tau,cr.nu);
    D = eta_scale(g,cr.P);
    rs = crn_stream(20260913,0,1); mdm = modem('4qam');
    X = gen_symbols(g.N,g.M,'4qam',rs); Xh = bsc_symbols(X,.08,mdm,rs);
    for ici = [false true]
        [A,~,Jt,b] = mcrb_AB_cp_ofdm(X,X,g,cr,ici,1,cr);
        Jg = fim(X,'cp-ofdm',g,cr,ici,1);
        pc = mcrb_assemble(A,Jt,b,Jg,cr.eta_bar,cr.eta_bar,g);
        lb = mcrb_eval(pc,.17);
        r = tc(r,sprintf('matched MCRB = LB = genie CRB, ICI %d',ici),max(relerr(lb.LB,lb.CRB),relerr(lb.MCRB0,lb.CRB)),1e-11);
        [A,B] = mcrb_AB_cp_otfs(X,Xh,g,cr,ici,.17,cr);
        [At,Bt] = mcrb_AB_cp_ofdm(isfft_dd(X),isfft_dd(Xh),g,cr,ici,.17,cr);
        r = tc(r,sprintf('DD/TF sandwich invariance, ICI %d',ici),relerr(stab_quad_bin(A,B,D),stab_quad_bin(At,Bt,D)),1e-12);
    end
    opts = struct('ftol',1e-12,'steptol',1e-10,'opttol',1e-10,'maxiter',500);
    [eta0,info] = pseudo_true(X,Xh,g,cr,true,opts);
    [a,t,n] = eta_unpack(eta0); ce = cr; ce.alpha=a; ce.tau=t; ce.nu=n;
    [A,~,Jt,b] = mcrb_AB_cp_ofdm(X,Xh,g,ce,true,1,cr);
    pc = mcrb_assemble(A,Jt,b,Jg,cr.eta_bar,eta0,g);
    r = tc(r,'pseudo-true projection converged',double(info.exitflag <= 0),0);
    f_ratio = norm(D\pc.F/D,'fro')/max(norm(D\pc.G/D,'fro'),eps);
    r = tc(r,'stationary pseudo-true score floor is negligible',f_ratio,1e-7);
    N0=.17;
    [An,Bn] = mcrb_AB_cp_ofdm(X,Xh,g,ce,true,N0,cr);
    lb = mcrb_eval(pc,N0);
    r = tc(r,'affine noise scaling agrees with direct sandwich',relerr(lb.MCRB0,stab_quad_bin(An,Bn,D)),1e-10);
    fprintf('    pseudo-true score-floor/noise ratio at N0=1: %.3e\n',f_ratio);
    opts.solver = 'lm';
    [eta_lm,info_lm] = pseudo_true(X,Xh,g,cr,true,opts);
    r = tc(r,'built-in LM converges to same projection residual',abs(info_lm.resnorm-info.resnorm)/info.resnorm,1e-9);
    r = tc(r,'built-in LM agrees in physical parameter bins',norm(D\(eta_lm-eta0)),1e-4);
end
