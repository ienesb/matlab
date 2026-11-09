function plot_hist_gain_constmag()
% PLOT_HIST_GAIN_CONSTMAG  Distribution of the per-subcarrier gain |H[n]| for the
% CONSTANT-MAGNITUDE reduced-EVA 3-tap channel (fixed tap amplitudes, random phase):
%   H[n] = sum_p amp_p * exp(-j2*pi*n*tap_p/N) * exp(j*phi_p),  phi_p ~ U(0,2pi).
% If |H[n]| ~ Rayleigh (and |H|^2 has finite density at 0 -> diversity order 1),
% that explains why OFDM's BER is order-1 (~1/SNR) even though the taps don't fade.
    thisdir=fileparts(mfilename('fullpath'));
    N=512; deltaf=15e3;
    eva_delays_ns=[0 150 370]; eva_pow_db=[0 -1.4 -0.6]; P=numel(eva_delays_ns);
    taps=round(eva_delays_ns*1e-9*N*deltaf);                 % [0 1 3]
    plin=10.^(eva_pow_db/20); amp=(plin/norm(plin)).';       % sum|amp|^2 = 1
    ns=(0:N-1).';
    ramp=exp(-1j*2*pi*ns*taps/N);                            % N x P subcarrier phase ramps

    rng(12345);
    nR=5000;                                                 % random-phase realizations
    Hall=zeros(N,nR);
    for r=1:nR
        Hall(:,r)=ramp*(amp.*exp(1j*2*pi*rand(P,1)));
    end
    H=Hall(:); g=abs(H); g2=g.^2;                            % ~2.5M samples
    fprintf('E|H|^2 = %.4f (should be 1);  E|H| = %.4f (Rayleigh: %.4f)\n', ...
            mean(g2), mean(g), sqrt(pi/4));

    hf=figure('Position',[160 120 980 460],'Color','w');

    % --- (a) pdf of |H| vs Rayleigh (E|H|^2 = 1) ---
    % subplot(1,2,1);
    histogram(g,'Normalization','pdf','NumBins',80,'FaceColor',[0.85 0.33 0.10], ...
              'EdgeColor','none','FaceAlpha',0.6); hold on; grid on; box on;
    rr=linspace(0,max(g),400);  fray=2*rr.*exp(-rr.^2);      % Rayleigh, E|H|^2=1
    plot(rr,fray,'k-','LineWidth',2);
    xlabel('|H[n]|'); ylabel('pdf');
    legend({'const-|amp| EVA (3 taps)','Rayleigh (E|H|^2=1)'},'Location','northeast','FontSize',9);
    title('Per-subcarrier gain |H[n]|  vs  Rayleigh');

    % --- (b) deep-fade tail: P(|H|^2 < eps) vs eps -> slope = diversity order ---
    % subplot(1,2,2);
    % eps_grid=logspace(-4,0,40);
    % Pfade=arrayfun(@(e) mean(g2<e), eps_grid);
    % loglog(eps_grid,Pfade,'o-','Color',[0.85 0.33 0.10],'LineWidth',1.8,'MarkerSize',4); hold on; grid on; box on;
    % loglog(eps_grid,eps_grid,'k--','LineWidth',1.5);          % slope 1 = order 1 (Rayleigh)
    % xlabel('\epsilon'); ylabel('P(|H|^2 < \epsilon)');
    % legend({'const-|amp| EVA','slope 1 (order-1 / Rayleigh)'},'Location','southeast','FontSize',9);
    % title('Deep-fade density  (slope = diversity order)');
    % 
    % sgtitle('Constant-magnitude reduced-EVA: per-subcarrier gain is Rayleigh-like (order 1)');
    if exist('theme','file'), theme(hf,"light"); end
    % figdir=fullfile(thisdir,'figures'); if ~exist(figdir,'dir'), mkdir(figdir); end
    % savefig(hf,fullfile(figdir,'hist_gain_constmag.fig'));
    % exportgraphics(hf,fullfile(thisdir,'results','hist_gain_constmag.png'),'Resolution',150,'BackgroundColor','white');
    % fprintf('Saved %s and .png\n',fullfile(figdir,'hist_gain_constmag.fig'));
    % 
    % % small-eps slope (empirical diversity order)
    % lo=eps_grid<1e-2 & Pfade>0;
    % pcoef=polyfit(log(eps_grid(lo)),log(Pfade(lo)),1);
    % fprintf('empirical small-eps slope (diversity order) = %.2f\n', pcoef(1));
end
