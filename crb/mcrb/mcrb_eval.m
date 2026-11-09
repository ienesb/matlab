function lb = mcrb_eval(pc, N0)
% MCRB_EVAL  Bounds at noise level N0 from the affine pieces of mcrb_assemble:
%   MCRB0 = N0 G + F,   LB = MCRB0 + BB,   CRB = N0 Gc
% plus per-block traces and RMSE bounds (sqrt of the trace) for the delay and
% Doppler blocks: lb.rmse_tau_{mcrb,lb,bb,crb}, lb.rmse_nu_*.
    idx = eta_index(pc.P);
    lb.MCRB0 = N0*pc.G + pc.F;
    lb.LB    = lb.MCRB0 + pc.BB;
    lb.CRB   = N0*pc.Gc;
    trd = @(X, i) sum(real(diag(X(i, i))));
    for blk = {'tau', 'nu'}
        i = idx.(blk{1});
        lb.(['tr_' blk{1} '_mcrb']) = trd(lb.MCRB0, i);
        lb.(['tr_' blk{1} '_bb'])   = trd(pc.BB, i);
        lb.(['tr_' blk{1} '_lb'])   = trd(lb.LB, i);
        lb.(['tr_' blk{1} '_crb'])  = trd(lb.CRB, i);
        lb.(['rmse_' blk{1} '_mcrb']) = sqrt(max(lb.(['tr_' blk{1} '_mcrb']), 0));
        lb.(['rmse_' blk{1} '_bb'])   = sqrt(max(lb.(['tr_' blk{1} '_bb']), 0));
        lb.(['rmse_' blk{1} '_lb'])   = sqrt(max(lb.(['tr_' blk{1} '_lb']), 0));
        lb.(['rmse_' blk{1} '_crb'])  = sqrt(max(lb.(['tr_' blk{1} '_crb']), 0));
    end
end
