function run_sim
% RUN_SIM  Reproduce the three CRLB figures of
%   B. Wang, J. Zhu, X. She, P. Chen, "Cramer-Rao Lower Bound Analysis for
%   OTFS and OFDM Modulation Systems," WPMC 2023.
%
% Four signal models (Sec. II): CP-Free OFDM, CP-OFDM, Zak-OTFS, Two-Steps OTFS.
% CRLB of delay (tau) and Doppler (nu) from the FIM of Sec. III, evaluated with
% a fixed random 4-QAM data block (nMC=1; raise NMC for Monte-Carlo averaging).
%
% Produces fig_1.png, fig_2.png, fig_3.png and results.mat in this folder.

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
cd(here);

df       = 15e3;                 % SCS = 15 kHz                     (Figs. 1-3)
% Channel matrices use the paper's literal per-scheme constants (Sec. III).
NMC      = 1;                    % fixed random realisation (upgrade -> e.g. 50)
SEED     = 2023;
schemes  = {'cpfree_ofdm','cp_ofdm','zak_otfs','twosteps_otfs'};
names    = {'CP-Free OFDM','CP-OFDM','Zak-OTFS','Two-Steps OTFS'};
col      = {[0.85 0.10 0.10],[0.95 0.55 0.00],[0.10 0.60 0.20],[0.10 0.30 0.90]};
sty      = {'-o','-s','-^','--x'};   % distinct per scheme (Zak/Two-Steps overlap)
SNR      = -20:4:20;             % dB                                (Figs. 1-3)

selftest(df);                   % sanity check on the channel matrices

% =====================================================================
% Figure 1 : single path, N=M=12, tau=3.33us, nu=500 Hz
% =====================================================================
M=12; N=12; tau=3.33e-6; nu=500;
sys = sys_params(M,N,df,tau);
path1 = struct('h',1,'tau',tau,'nu',nu);
F1 = struct();
i0 = find(SNR==0,1);
for s = 1:numel(schemes)
    r = compute_crlb_scheme(schemes{s}, path1, sys, SNR, NMC, SEED);
    F1.(schemes{s}) = r;
    fprintf('Fig1 %-15s cond(Gt)=%.2e  CRLB@0dB: tau=%.3e s^2  nu=%.3e Hz^2\n', ...
        names{s}, r.condGt, r.tau(1,i0), r.nu(1,i0));
end
plot_two(SNR, F1, schemes, names, col, sty, ...
    'Fig. 1: CRLB vs SNR  (N=M=12, \tau=3.33\mus, \nu=500 Hz)');
saveas(gcf, fullfile(here,'fig_1.png'));

% =====================================================================
% Figure 2 : CRLB vs SNR for several sizes M=N (other params as Fig. 1)
% =====================================================================
sizes = [8 12 16];
ls    = {':','-','--'};                      % linestyle per size
figure('Position',[100 100 1100 460],'Color','w');
for sp = 1:2                                  % 1=delay, 2=Doppler
    ax = subplot(1,2,sp); hold on; grid on; box on; set(ax,'Color','w');
    for s = 1:numel(schemes)
        for q = 1:numel(sizes)
            Mq = sizes(q); Nq = sizes(q);
            sysq = sys_params(Mq,Nq,df,tau);
            r = compute_crlb_scheme(schemes{s}, ...
                    struct('h',1,'tau',tau,'nu',nu), sysq, SNR, NMC, SEED);
            y = (sp==1).*r.tau + (sp==2).*r.nu;
            semilogy(SNR, y, ls{q}, 'Color', col{s}, 'LineWidth',1.4);
        end
    end
    set(ax,'YScale','log'); xlabel('SNR [dB]');
    if sp==1, ylabel('CRLB of \tau  [s^2]'); title('(a) delay vs SNR');
    else,     ylabel('CRLB of \nu  [Hz^2]'); title('(b) Doppler vs SNR'); end
end
F2.sizes=sizes; F2.ls=ls;
legend_sizes(names, col, sizes, ls);
sgtitle('Fig. 2: CRLB vs SNR for sizes M=N \in \{8,12,16\}');
saveas(gcf, fullfile(here,'fig_2.png'));

% =====================================================================
% Figure 3 : two-path channel, N=M=12
%   h1=0.7e^{j pi/3}, h2=0.3e^{j 3pi/4}
%   tau1=3.33us nu1=500Hz ,  tau2=5us nu2=2.5kHz
% =====================================================================
M=12; N=12;
paths(1) = struct('h',0.7*exp(1j*pi/3),  'tau',3.33e-6,'nu',500);
paths(2) = struct('h',0.3*exp(1j*3*pi/4),'tau',5e-6,   'nu',2500);
sys = sys_params(M,N,df,max([paths.tau]));
F3 = struct();
for s = 1:numel(schemes)
    r = compute_crlb_scheme(schemes{s}, paths, sys, SNR, NMC, SEED);
    F3.(schemes{s}) = r;
    fprintf('Fig3 %-15s cond(Gt)=%.2e  CRLB@0dB: t1=%.2e t2=%.2e  n1=%.2e n2=%.2e\n', ...
        names{s}, r.condGt, r.tau(1,i0), r.tau(2,i0), r.nu(1,i0), r.nu(2,i0));
end
plot_two_2path(SNR, F3, schemes, names, col, ...
    'Fig. 3: CRLB vs SNR, two-path channel (path 1: \tau=3.33\mus,\nu=500Hz; path 2: \tau=5\mus,\nu=2.5kHz)');
saveas(gcf, fullfile(here,'fig_3.png'));

save(fullfile(here,'results.mat'),'SNR','F1','F2','F3','schemes','names','df','NMC');
fprintf('\nSaved fig_1.png, fig_2.png, fig_3.png and results.mat\n');
end

% =====================================================================
function plot_two(SNR, F, schemes, names, col, sty, ttl)
figure('Position',[100 100 1100 460],'Color','w');
for sp = 1:2
    ax = subplot(1,2,sp); hold on; grid on; box on; set(ax,'Color','w');
    for s = 1:numel(schemes)
        r = F.(schemes{s});
        y = (sp==1).*r.tau(1,:) + (sp==2).*r.nu(1,:);
        semilogy(SNR, y, sty{s}, 'Color',col{s},'LineWidth',1.6,'MarkerSize',5);
    end
    set(ax,'YScale','log'); xlabel('SNR [dB]');
    if sp==1, ylabel('CRLB of \tau  [s^2]'); title('(a) delay vs SNR');
    else,     ylabel('CRLB of \nu  [Hz^2]'); title('(b) Doppler vs SNR'); end
    if sp==1, legend(names,'Location','southwest','Color','w','TextColor','k'); end
end
sgtitle(ttl);
end

% =====================================================================
function plot_two_2path(SNR, F, schemes, names, col, ttl)
figure('Position',[100 100 1100 460],'Color','w');
mk = {'-o','--s'};   % path 1 solid, path 2 dashed
for sp = 1:2
    ax = subplot(1,2,sp); hold on; grid on; box on; set(ax,'Color','w');
    leg = {};
    for s = 1:numel(schemes)
        r = F.(schemes{s});
        for p = 1:2
            y = (sp==1).*r.tau(p,:) + (sp==2).*r.nu(p,:);
            semilogy(SNR, y, mk{p}, 'Color',col{s},'LineWidth',1.4,'MarkerSize',4);
            leg{end+1} = sprintf('%s, path %d', names{s}, p); %#ok<AGROW>
        end
    end
    set(ax,'YScale','log'); xlabel('SNR [dB]');
    if sp==1, ylabel('CRLB of \tau  [s^2]'); title('(a) delay vs SNR');
    else,     ylabel('CRLB of \nu  [Hz^2]'); title('(b) Doppler vs SNR'); end
    if sp==1, legend(leg,'Location','southwest','FontSize',7,'NumColumns',2,'Color','w','TextColor','k'); end
end
sgtitle(ttl);
end

% =====================================================================
function legend_sizes(names, col, sizes, ls)
% Build a compact two-part legend (scheme colours + size linestyles).
h = []; lab = {};
for s = 1:numel(names)
    h(end+1) = plot(nan,nan,'-','Color',col{s},'LineWidth',2); %#ok<AGROW>
    lab{end+1} = names{s}; %#ok<AGROW>
end
for q = 1:numel(sizes)
    h(end+1) = plot(nan,nan,ls{q},'Color',[.4 .4 .4],'LineWidth',1.5); %#ok<AGROW>
    lab{end+1} = sprintf('M=N=%d',sizes(q)); %#ok<AGROW>
end
legend(h,lab,'Location','southwest','FontSize',7,'Color','w','TextColor','k');
end

% =====================================================================
function selftest(df)
% Report the Frobenius energy of each scheme's closed-form Psi for a single
% on-grid path (h=1, tau=3.33us, nu=500Hz, N=M=12).  Sanity check that the
% paper's per-scheme global constants are in force and matrices are finite.
M=12; N=12; tau=3.33e-6; nu=500;
sys = sys_params(M,N,df,tau);
for sc = {'cpfree_ofdm','cp_ofdm','zak_otfs','twosteps_otfs'}
    Psi = psi_build(sc{1}, 1, tau, nu, sys);
    fprintf('selftest %-15s  ||Psi||_F^2 = %.3e\n', sc{1}, norm(Psi,'fro')^2);
end
end
