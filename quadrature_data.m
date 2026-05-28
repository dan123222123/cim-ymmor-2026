%% Setup
n = 6; refew = -n:-1;
A = diag(refew); B = ones(n,2); C = ones(n,2)';
%
H = @(z) C*((z*eye(size(A)) - A) \ B);
%
% Virginia Tech palette (matches tex/main.tex preamble: RGB(134,31,65) and
% RGB(232,119,34)). MATLAB's scatter() rejects hex strings as positional
% color, so keep these as RGB triplets.
VTmaroon = [134  31  65]/255;
VTorange = [232 119  34]/255;
%
% Figure-save directory: ../tex/figures/ relative to this script. Resolve
% once so cells run later still know where to write.
thisfile = mfilename('fullpath');
if isempty(thisfile); thisfile = fullfile(pwd,'quadrature_data.m'); end
figDir = fullfile(fileparts(thisfile),'..','tex','figures');
if ~exist(figDir,'dir'); mkdir(figDir); end

%% Quadrature Data, MIMO Case
import Visual.*; % allows us to skip subsequent "Visual."s

% say we take an ellipse about our spectrum
% to do this, we specify the center, horizontal/vertical semiradii, and the
% initial number of quadrature nodes; the N below is overwritten per pass in
% the loop further down.
c = Contour.Ellipse(-(n+1)/2,n/3,n/4,16); c.plot_quadrature = true;
% specify our operator of interest
o = OperatorData(H);
% along with its poles and sampling mode
o.refew = diag(A); o.sample_mode = "Direct";

% In-contour reference eigenvalues, hoisted once and forced to column. The
% contour geometry doesn't change inside the per-N sweep (only c.N does), so
% we can compute the inside-set once and broadcast a plain vector into the
% parfor -- avoids calling c.inside (a handle method) on every worker, and
% keeps the shape compatible with eig()'s column-vector output below.
inside_ew = refew(c.inside(refew)).';
% Order of the local realization: rank-m SVD truncation of the K x K data
% matrix Db isolates the m eigenvalues inside the contour. See the slide on
% "Practical Notes: Quadrature and Oversampling" for the math.
m = length(inside_ew);

% with the operator and sampling data set, we can initialize our sampling
% data structure, and assign a plotting axis associated to the operator
% reference poles and contour.
figure(1); clf;
s = SampleData(o,c); s.ax = gca;
axis equal; grid on;
xlabel('$\mathrm{Re}\,z$','Interpreter','latex');
ylabel('$\mathrm{Im}\,z$','Interpreter','latex');
title('Ellipse contour with quadrature nodes','Interpreter','latex');

% sketch evaluations of our operator down to the SISO case so the moment
% builds stay light during the per-N grid sweep.
s.ell = 1; s.r = 1;

%% Set Data Matrix Size & heatmap grid (shared across all N)
K = n;
Ngrid = 1001; xg = linspace(-10,5,Ngrid); yg = linspace(-7.5,7.5,Ngrid);
[Xg,Yg] = meshgrid(xg,yg); G = Xg + 1i*Yg;

%% Sweep over a few quadrature-node counts.
N_values = [8, 32, 256];
nN = numel(N_values);
%
% K rows: SVD of the full K x K data matrix, so the decay plot exposes the
% rank-m gap (sigma_{m+1}/sigma_1 should crash) rather than truncating to m
% values up-front.
Sigma_era_all   = zeros(K, nN);
Sigma_spl_all   = zeros(K, nN);
sigma_best_all  = zeros(1, nN);
ERA_err_all     = zeros(1, nN);
best_spl_err_all= zeros(1, nN);
% Cache the per-N heatmap arrays + quadrature nodes so the rendering pass
% below can build a single 3-panel montage with a shared colormap/colorbar.
ls_eravspl_all  = cell(1, nN);
z_all           = cell(1, nN);
%
for idx = 1:nN
    Nq = N_values(idx);
    % refresh quadrature at this node count
    c.N = Nq; s.compute();
    z = c.z; w = c.w; Ql = s.Ql; Qr = s.Qr; Qlr = s.Qlr;

    % ERA (sigma = Inf). build_sploewner now returns the moment block B/C
    % too -- we feed those plus Db, Ds to Numerics.realize, which does the
    % rank-m SVD truncation: Db = X Sigma Y', keep top m columns, then solve
    % the m-by-m GEP eig(sqrt(Sigma)\(X'*Ds*Y)/sqrt(Sigma)). Returns exactly
    % m eigenvalues -- the in-contour ones -- so gmatch pairs 1-1 with
    % inside_ew. Sigma_full (= Sigmaf) is the full K-by-K SVD diagonal we
    % want for the decay plot.
    sigma = Inf;
    [Ml,Mr,Mlr] = Numerics.sploewner.build_quadrature_moments(sigma,z,w,Ql,Qr,Qlr,K);
    [Db,Ds,B,C] = Numerics.sploewner.build_sploewner(sigma,Ml,Mr,Mlr,K);
    try
        [ew,~,~,~,Sigma_full] = Numerics.realize(m,Db,Ds,B,C,NaN);
        ERA_err = gmatch(inside_ew, ew);
        Sigma_era = diag(Sigma_full); Sigma_era = Sigma_era / Sigma_era(1);
    catch e
        fprintf(2,'Warning: ERA realization failed at N=%d: %s\n',Nq,e.message);
        ERA_err = NaN;
        Sigma_era = nan(K,1);
    end

    % SPLoewner over the sigma grid. Numerics.realize errors if the K x K
    % data matrix has numerical rank < m (which happens when sigma sits on or
    % very near a pole of the realized pencil); catch and leave NaN so the
    % heatmap renders those points as transparent holes.
    SPLoewner_err = nan(Ngrid,Ngrid);
    parfor i=1:Ngrid
        for j=1:Ngrid
            sigma_loc = G(i,j);
            try
                [Ml_l,Mr_l,Mlr_l] = Numerics.sploewner.build_quadrature_moments( ...
                                           sigma_loc,z,w,Ql,Qr,Qlr,K);
                [Db_l,Ds_l,B_l,C_l] = Numerics.sploewner.build_sploewner( ...
                                           sigma_loc,Ml_l,Mr_l,Mlr_l,K);
                ew_l = Numerics.realize(m,Db_l,Ds_l,B_l,C_l,NaN);
                if all(isfinite(ew_l))
                    SPLoewner_err(i,j) = gmatch(inside_ew, ew_l);
                end
            catch
                % leave as NaN
            end
        end
    end

    ls_eravspl = log10(SPLoewner_err./ERA_err);
    [~,bsidx] = min(ls_eravspl,[],"all"); sigma_best = G(bsidx);

    % SPLoewner SVD at the best sigma for this N. Same realize path as ERA;
    % we only need the K-vector of singular values for the decay plot.
    [Ml,Mr,Mlr] = Numerics.sploewner.build_quadrature_moments(sigma_best,z,w,Ql,Qr,Qlr,K);
    [Db,Ds,B,C] = Numerics.sploewner.build_sploewner(sigma_best,Ml,Mr,Mlr,K);
    try
        [~,~,~,~,Sigma_full] = Numerics.realize(m,Db,Ds,B,C,NaN);
        Sigma_spl = diag(Sigma_full); Sigma_spl = Sigma_spl / Sigma_spl(1);
    catch
        Sigma_spl = nan(K,1);
    end

    % Stash for the montage render + joint .dat / meta.tex writeout below.
    ls_eravspl_all{idx}     = ls_eravspl;
    z_all{idx}              = z;
    Sigma_era_all(:,idx)    = Sigma_era;
    Sigma_spl_all(:,idx)    = Sigma_spl;
    sigma_best_all(idx)     = sigma_best;
    ERA_err_all(idx)        = ERA_err;
    best_spl_err_all(idx)   = SPLoewner_err(bsidx);
end

%% 3-panel heatmap montage with a shared colorbar.
% The slide ("Quadrature Error Changes the Picture" in tex/main.tex) loads a
% single PNG so each panel inherits the same color scale + a single legible
% colorbar; per-panel legends are intentionally omitted (marker meaning is
% covered by the preceding "Practical Notes" slide).
figure(2); clf;
% Portrait aspect: 3 square data panels stacked + ~10% width colorbar.
set(gcf, 'Position', [100 100 540 1280]);
tl = tiledlayout(3,1,'Padding','compact','TileSpacing','compact');
for idx = 1:nN
    Nq = N_values(idx);
    nexttile;
    % AlphaData masks NaN cells transparent so the axis Color punches through
    % -- otherwise NaN renders against the default white background and gets
    % visually confused with mid-colormap values near log10(SPL/ERA) = 0.
    A = ls_eravspl_all{idx};
    him = imagesc(xg,yg,A,'AlphaData',~isnan(A));
    axis xy; axis equal;
    set(gca,'Color','g');                  % NaN cells -> opaque green
    xlim([xg(1) xg(end)]); ylim([yg(1) yg(end)]);
    clim([-1 1]);
    hold on;
    % Ellipse boundary
    t = linspace(0,2*pi,512);
    ellipse_z = c.gamma + c.alpha*cos(t) + 1i*c.beta*sin(t);
    plot(real(ellipse_z),imag(ellipse_z),'-','Color',VTmaroon,'LineWidth',1.5);
    % Quadrature nodes (only when sparse enough to read)
    if Nq < 64
        scatter(real(z_all{idx}),imag(z_all{idx}), 40, VTmaroon, 'x','LineWidth',1.2);
    end
    % Reference eigenvalues + best-sigma marker
    scatter(real(refew),imag(refew), 70, VTorange, 'd', 'filled', ...
            'MarkerEdgeColor','k','LineWidth',0.8);
    scatter(real(sigma_best_all(idx)),imag(sigma_best_all(idx)), 90, 'k', 's', 'filled', ...
            'MarkerEdgeColor','w','LineWidth',0.8);
    hold off;
    title(sprintf('$N=%d$\nERA err $=%.2e$\nbest SPL err $=%.2e$', ...
                  Nq, ERA_err_all(idx), best_spl_err_all(idx)), ...
          'Interpreter','latex','FontSize',11);
    ylabel('$\mathrm{Im}\,\sigma$','Interpreter','latex');
    if idx == nN
        xlabel('$\mathrm{Re}\,\sigma$','Interpreter','latex');
    else
        set(gca,'XTickLabel',[]);
    end
    set(gca,'Layer','top');
end
colormap(redblue(5000));
cb = colorbar; cb.Layout.Tile = 'east';
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\log_{10}(\|\mathrm{SPLoewner}\| / \|\mathrm{ERA}\|)$';
exportgraphics(gcf, fullfile(figDir,'quad_heatmaps.png'), 'Resolution', 500);

%% Singular-value decay -- exported for pgfplots.
% The "Quadrature Error Changes the Picture" frame in tex/main.tex renders
% these curves in LaTeX so the math fonts match the rest of the deck. Two
% files written:
%   tex/figures/quad_svd.dat       -- one row per k, columns
%                                     k era<N1> spl<N1> ... era<Nk> spl<Nk>
%                                     where <Ni> are the entries of N_values
%   tex/figures/quad_svd_meta.tex  -- defines \QuadSvdSplLabelN<N> via
%                                     \csname so the legend tracks sigma_best
% The MATLAB plotting code below is intentionally left commented so a user can
% uncomment it and preview the curves locally.
%
% Build a single multi-column table; pgfplots' "table[x=k, y=era<N>]" picks
% columns by header name.
varNames = cell(1, 1 + 2*nN);
varNames{1} = 'k';
data = cell(1, 1 + 2*nN);
% Sigma_era_all is (K x nN), not (n x nN), since K may differ from n when we
% only target the in-contour spectrum -- match the table's row count to the
% actual SVD length so writetable() doesn't reject mismatched columns.
data{1} = (1:size(Sigma_era_all,1)).';
for idx = 1:nN
    Nq = N_values(idx);
    varNames{2*idx}   = sprintf('era%d', Nq);
    varNames{2*idx+1} = sprintf('spl%d', Nq);
    data{2*idx}       = Sigma_era_all(:,idx);
    data{2*idx+1}     = Sigma_spl_all(:,idx);
end
T = table(data{:}, 'VariableNames', varNames);
writetable(T, fullfile(figDir,'quad_svd.dat'), ...
           'Delimiter','space', 'FileType','text');
%
% Write the legend-label macros. Macro names contain digits, so we go through
% \csname; consumer side: \csname QuadSvdSplLabelN16\endcsname.
fid = fopen(fullfile(figDir,'quad_svd_meta.tex'),'w');
fprintf(fid, '%% Auto-generated by code/quadrature_data.m -- do not edit by hand.\n');
for idx = 1:nN
    Nq = N_values(idx); sb = sigma_best_all(idx);
    fprintf(fid, ['\\expandafter\\providecommand\\csname QuadSvdSplLabelN%d\\endcsname' ...
                  '{SPLoewner ($N=%d$, $\\sigma=%.2f%+0.2fi$)}\n'], ...
            Nq, Nq, real(sb), imag(sb));
end
fclose(fid);

%% Local MATLAB preview of the SVD decay (commented out; the slide uses pgfplots
%% via the .dat/.tex files written above. Uncomment this block to preview in MATLAB).
%{
figure(2+nN); clf;
markers = {'-o','-s','-^'};
kk = (1:size(Sigma_era_all,1)).';
hold on;
for idx = 1:nN
    Nq = N_values(idx); sb = sigma_best_all(idx);
    plot(kk,Sigma_era_all(:,idx),markers{idx},'LineWidth',1.5, ...
         'MarkerFaceColor',VTmaroon,'Color',VTmaroon, ...
         'DisplayName',sprintf('ERA ($N=%d$)',Nq));
    plot(kk,Sigma_spl_all(:,idx),markers{idx},'LineWidth',1.5, ...
         'MarkerFaceColor',VTorange,'Color',VTorange, ...
         'DisplayName',sprintf('SPLoewner ($N=%d$, $\\sigma=%.2f%+0.2f\\,i$)', ...
                               Nq,real(sb),imag(sb)));
end
hold off;
yscale("log"); xlim([kk(1),kk(end)]); xticks(kk); grid on;
xlabel('$k$','Interpreter','latex');
ylabel('$\sigma_k/\sigma_1$','Interpreter','latex');
legend('Interpreter','latex','Location','northoutside','Orientation','vertical','FontSize',8);
exportgraphics(gcf, fullfile(figDir,'quad_svd.png'), 'Resolution', 500);
%}
%
for idx = 1:nN
    fprintf("N=%d -- ERA Error: %e vs Best SPLoewner Error %e\n", ...
            N_values(idx), ERA_err_all(idx), best_spl_err_all(idx));
end
