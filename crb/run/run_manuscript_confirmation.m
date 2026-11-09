function run_manuscript_confirmation()
% Fresh-seed confirmation at the two low-SNR comparison points. This writes
% a new result family and refuses to overwrite either completed or partial
% outputs. Run from matlab/crb after setup_crb.
    root = fileparts(fileparts(mfilename('fullpath')));
    tag = '_finalization20260914_seed9101';
    variants = {'reduced','reduced-fast'};
    for k = 1:numel(variants)
        channel = ['eva_' strrep(variants{k},'-','_')];
        names = {fullfile(root,'results',['mcrb_rx_' channel '_n100' tag '.mat']), ...
            fullfile(root,'results',['CKPT_mcrb_rx_' channel '_n100' tag '.mat']), ...
            fullfile(root,'figures',['fig_mcrb_rx_' channel '_n100' tag '.png']), ...
            fullfile(root,'figures',['fig_ber_mcrb_panels_' channel '_n100' tag '.png'])};
        for j = 1:numel(names)
            assert(~isfile(names{j}),'run_manuscript_confirmation:exists', ...
                'Refusing to overwrite existing output: %s',names{j});
        end
    end
    run_manuscript_checks();
    if isempty(gcp('nocreate')), parpool('Processes',8); end
    t0 = tic;
    for k = 1:numel(variants)
        fprintf('\nCONFIRMATION variant=%s, seed=9101, trials=100, SNR=[4 8]\n',variants{k});
        r = run_mcrb_rx_eva('tier',100,'variant',variants{k},'snr_db',[4 8], ...
            'seed',9101,'tag',tag,'bsc_mode','uniform');
        fprintf('CONFIRMATION_COMPLETED variant=%s seconds=%.3f elapsed=%.3f\n', ...
            variants{k},r.seconds,toc(t0));
    end
    fprintf('ALL_CONFIRMATIONS_COMPLETED elapsed=%.3f\n',toc(t0));
end
