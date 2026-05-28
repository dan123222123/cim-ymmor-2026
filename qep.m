import Visual.*;

% From Tisseur&Meerbergen 2001, section 3.6

M = [1 0 0;0 1 0;0 0 0];
C = [-2 0 1;0 0 0;0 0 0];
K = [1 0 0;0 -1 0;0 0 1];
T = @(z) z.^2*M + z.*C + K;

% Virginia Tech palette (matches tex/main.tex preamble: RGB(134,31,65) and
% RGB(232,119,34)). MATLAB's scatter() rejects hex strings as positional
% color, so keep these as RGB triplets; StylePreferences wants the hex form.
VTmaroon = [134  31  65]/255;
VTorange = [232 119  34]/255;
VTmaroon_hex = '#861F41';
VTorange_hex = '#E87722';

% Figure-save directory.
thisfile = mfilename('fullpath');
if isempty(thisfile); thisfile = fullfile(pwd,'qep.m'); end
figDir = fullfile(fileparts(thisfile),'..','tex','figures');
if ~exist(figDir,'dir'); mkdir(figDir); end

n = OperatorData(T,'qep_3.6');
c = Contour.Ellipse(0,2,0.5);
cim = CIM(n,c);

cim.SampleData.OperatorData.refew = [-1;1;1;1;Inf;Inf];

% Style-preferences theme matching the deck (Virginia Tech maroon + orange).
% This same struct interface is documented in CIMTOOL/src/demos/style_preferences_demo.m.
vtTheme = struct( ...
    'ContourColor',                VTmaroon_hex, ...
    'ContourLineWidth',            4, ...
    'QuadratureColor',             VTmaroon_hex, ...
    'QuadratureMarker',            'x', ...
    'ShowQuadratureNodes',         true, ...
    'ReferenceEigenvalueColor',    VTorange_hex, ...
    'ReferenceEigenvalueMarker',   'diamond', ...
    'ComputedEigenvalueColor',     '#FFD700', ...   % gold, against maroon contour
    'ComputedEigenvalueMarker',    'o', ...
    'SPLoewnerShiftColor',         'k', ...
    'SPLoewnerShiftMarker',        'square' ...
);

cimtool = CIMTOOL(cim, vtTheme);

%% some text modes
cim.setComputationalMode(Numerics.ComputationalMode.Hankel);

% cim.SampleData.ell = 2; cim.SampleData.r = 2;
% cim.RealizationData.RealizationSize = Numerics.RealizationSize(4,2);
% cim.SampleData.Contour.N = 128;

cim.SampleData.ell = 2; cim.SampleData.r = 2;
cim.RealizationData.RealizationSize = Numerics.RealizationSize(4,8);
cim.SampleData.Contour.N = 64;

ew = cim.eigs(); ERA_err = gmatch([-1;1;1;1],ew)

% %% ERA vs SPLoewner Heatmap
% cim.setComputationalMode(Numerics.ComputationalMode.SPLoewner); cim.eigs;
% % The parfor below copies a handle object per worker. Two requirements:
% %   (1) use a *process* pool -- thread pools (the R2025b default) cannot
% %       copy() handle classes, so copy(ncim) errors on a thread worker;
% %   (2) broadcast a graphics-free Numerics.CIM via toNumerics() -- shipping
% %       the Visual.CIM (scatter handles + property listeners) otherwise spews
% %       "proplistener"/"Scatter" load warnings on every worker.
% pp = gcp('nocreate');
% if isempty(pp) || ~isa(pp,'parallel.ProcessPool')
%     if ~isempty(pp); delete(pp); end
%     parpool('Processes');
% end
% ncim = cim.toNumerics();
% N = 51; x = linspace(-5,5,N); [X,Y] = meshgrid(x,x); G = X + 1i*Y;
% SPLoewner_err = zeros(N,N);
% parfor i=1:N
%     ccim = copy(ncim);
%     for j=1:N
%         sigma = G(i,j);
%         ccim.RealizationData.InterpolationData = Numerics.InterpolationData([],sigma);
%         try
%             ew = ccim.eigs; SPLoewner_err(i,j) = gmatch([-1;1;1;1],ew);
%         catch e
%             SPLoewner_err(i,j) = NaN;
%         end
%     end
% end
% %% Plot heatmap with VT-themed contour/refew/best-sigma overlay
% ls_eravspl = log10(SPLoewner_err./ERA_err);
% [bsn,bsidx] = min(ls_eravspl,[],"all"); sigma_best = G(bsidx);
% refew_qep = [-1;1;1;1];
% %
% figure(2); clf;
% imagesc(x,x,ls_eravspl); axis xy; axis equal;
% xlim([x(1) x(end)]); ylim([x(1) x(end)]);
% colormap(redblue(5000)); clim([-1 1]);
% cb = colorbar; cb.Label.Interpreter = 'latex';
% cb.Label.String = '$\log_{10}(\|\mathrm{SPLoewner}\| / \|\mathrm{ERA}\|)$';
% xlabel('$\mathrm{Re}\,\sigma$','Interpreter','latex');
% ylabel('$\mathrm{Im}\,\sigma$','Interpreter','latex');
% title(sprintf('Tisseur \\& Meerbergen \\S3.6 QEP: ERA err $= %.2e$;\\ best SPLoewner err $= %.2e$', ...
%               ERA_err, SPLoewner_err(bsidx)), 'Interpreter','latex');
% hold on;
% % Ellipse boundary
% ec = cim.SampleData.Contour;
% t = linspace(0,2*pi,512);
% ellipse_z = ec.gamma + ec.alpha*cos(t) + 1i*ec.beta*sin(t);
% plot(real(ellipse_z),imag(ellipse_z),'-','Color',VTmaroon,'LineWidth',1.5, ...
%      'DisplayName','$\partial\mathcal{D}$');
% % Reference eigenvalues (the four finite ones in [-1, 1])
% scatter(real(refew_qep),imag(refew_qep), 100, VTorange, 'd', 'filled', ...
%         'MarkerEdgeColor','k','LineWidth',1.0, ...
%         'DisplayName','exact eigenvalues');
% % Best-sigma marker
% scatter(real(sigma_best),imag(sigma_best), 120, 'k', 's', 'filled', ...
%         'MarkerEdgeColor','w','LineWidth',1.0, ...
%         'DisplayName','best $\sigma$');
% hold off;
% legend('Interpreter','latex','Location','northeast','Color','w','FontSize',8);
% set(gca,'Layer','top');
% exportgraphics(gcf, fullfile(figDir,'qep_heatmap.png'), 'Resolution', 220);
% 
% %% SVD decay at sigma=Inf vs best sigma
% figure(3); clf; tl = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
% nexttile;
% cim.setComputationalMode(Numerics.ComputationalMode.Hankel);
% cim.eigs; Sigma = diag(cim.ResultData.Sigma); Sigma = Sigma / Sigma(1);
% plot(1:length(Sigma),Sigma,'-o','LineWidth',1.5, ...
%      'MarkerFaceColor',VTmaroon,'Color',VTmaroon);
% yscale("log"); xlim([1,length(Sigma)]); grid on;
% xlabel('$k$','Interpreter','latex');
% ylabel('$\sigma_k/\sigma_1$','Interpreter','latex');
% title('ERA ($\sigma=\infty$)','Interpreter','latex');
% 
% nexttile;
% cim.setComputationalMode(Numerics.ComputationalMode.SPLoewner);
% cim.RealizationData.InterpolationData = Numerics.InterpolationData([],sigma_best);
% cim.eigs; Sigma = diag(cim.ResultData.Sigma); Sigma = Sigma / Sigma(1);
% plot(1:length(Sigma),Sigma,'-o','LineWidth',1.5, ...
%      'MarkerFaceColor',VTmaroon,'Color',VTmaroon);
% yscale("log"); xlim([1,length(Sigma)]); grid on;
% xlabel('$k$','Interpreter','latex');
% title(sprintf('SPLoewner at best $\\sigma=%.2f%+0.2f\\,i$', ...
%               real(sigma_best),imag(sigma_best)), 'Interpreter','latex');
% title(tl,'Normalized singular-value decay of $\mathbf{D}$', 'Interpreter','latex');
% exportgraphics(gcf, fullfile(figDir,'qep_svd.png'), 'Resolution', 220);
% %
% fprintf("ERA Error: %e vs Best SPLoewner Error %e\n",ERA_err,SPLoewner_err(bsidx))
