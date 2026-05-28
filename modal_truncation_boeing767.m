% MODAL_TRUNCATION_BOEING767  Demo 4b: stable/unstable split, Boeing 767 SISO.
%
%   The Boeing 767 aeroelastic model (n=55) has a mix of stable and unstable
%   poles -- a handful of right-half-plane modes drive the unstable behaviour
%   of the elevator-to-pitch-rate channel. CIMTOOL's Numerics.ModalTruncation
%   isolates the unstable part directly from frequency-domain samples of the
%   SISO transfer function H(z) = c^T (zI - A)^{-1} b, returning callable
%   handles for H_D (the in-contour unstable subsystem) and the residual
%   H_{D^c} = H - H_D.
%
%   Produces two figures referenced by Demo 4b of the YMMOR deck:
%       tex/figures/mt_boeing_poles.pdf -- contour + true vs. computed poles
%       tex/figures/mt_boeing_bode.pdf  -- |H|, |H_D|, |H_{D^c}| vs. omega
%
%   The script runs in two modes (see "Run mode" cell below):
%       run_mode = "demo":  fast defaults -- [64;64] quadrature nodes, suitable
%                           for a live cell-by-cell demo in MATLAB.
%       run_mode = "figs":  larger N for the converged figures used on the slide.
%
%   Requires CIMTOOL (math/CIMTOOL) on the MATLAB path; the Paths cell below
%   adds it relative to this file. The Boeing 767 .mat ships alongside this
%   script in code/data/. We define H from the original (A,B,C) directly --
%   the diagonalising similarity transform used elsewhere has condition number
%   ~1e20 and is not safe to sample through.
%
%   See also: Numerics.ModalTruncation, Numerics.Contour.CircularSegment,
%             code/modal_truncation.m (synthetic counterpart, Demo 4a)

%% Paths
here = fileparts(mfilename('fullpath'));
if exist('Numerics.ModalTruncation','class') ~= 8
    addpath(genpath(fullfile(here,'..','..','..','..','CIMTOOL','src')));
end
% Data ships in code/data/ so this folder is self-contained for a future repo.
data_dir = fullfile(here,'data');
fig_dir  = fullfile(here,'..','tex','figures');
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
fsp = 16;  % figure font size

%% Run mode: "demo" (live) or "figs" (slide artifacts)
% Cell-by-cell users typically run the script as-is in "demo" mode and only
% flip to "figs" when the slide figures need regenerating. Pre-setting
% `run_mode` in the workspace before calling this script (e.g. from a
% batch wrapper) is also honoured.
if ~exist('run_mode','var') || isempty(run_mode)
    run_mode = "demo";
end
switch run_mode
    case "demo"
        Nquad = [64; 64];      % nodes on [arc; chord]
    case "figs"
        Nquad = [256; 256];    % converged contour for the slide PDF
    otherwise
        error('Unknown run_mode "%s" (use "demo" or "figs").', run_mode);
end

%% Load Boeing 767 SISO model (ELEV -> QCG channel)
% The .mat ships A,B,C,n,m,p for the full MIMO model; we slice to SISO on the
% elevator-to-pitch-rate channel, matching ddmt/models/boeing767_nnLTI_SISO.m.
S = load(fullfile(data_dir,'boeing767_nnLTI.mat'),'A','B','C');
A = S.A; b = S.B(:,1); c = S.C(1,:);
n = size(A,1);

% Reference eigenvalues -- partition by stability for the Figure 1 scatter.
ewref     = eig(A);
ew_ust    = ewref(real(ewref) >  0);
ew_stab   = ewref(real(ewref) <= 0);
fprintf('Boeing 767 SISO: n = %d, #unstable = %d, #stable = %d\n', ...
        n, numel(ew_ust), numel(ew_stab));

% SISO transfer function as a callable. Defining H from (A,b,c) directly avoids
% the rank-deficient diagonalising similarity used in pole-residue form (the
% diagonaliser V has cond(V) ~ 1e20 for this model).
I_n = eye(n);
H   = @(z) c * ((z*I_n - A) \ b);

%% Circular-segment contour over {Re z > 0}, enclosing the unstable poles
% A circle of radius rho centered at gamma=0 covering the right half-plane
% (theta = [-pi/2, pi/2]) bounds {Re z > 0} ∩ {|z| < rho}. Sized so the small
% cluster of unstable Boeing poles sits comfortably inside.
rho_ust = max(abs(ew_ust));         % furthest unstable pole from the origin
rho     = 1.5 * rho_ust;            % a generous bracket around the cluster
gamma   = 0;
phi     = [-pi/2, pi/2];
contour = Numerics.Contour.CircularSegment(gamma, rho, phi, Nquad);

%% Modal truncation: realize the in-region (unstable) subsystem via CIMTOOL
% Numerics.ModalTruncation wires up an OperatorData/SampleData/CIM/Realization
% pipeline behind the scenes; we just feed it H, the contour, and the size of
% the in-region system we expect (#unstable poles).
m_in  = numel(ew_ust);              % expected #poles inside D
T_blk = max(m_in + 1, 5);           % block-Hankel rows/cols (sketching dim)

rd = Numerics.RealizationData();
rd.RealizationSize   = Numerics.RealizationSize(m_in, T_blk, T_blk);
rd.ComputationalMode = Numerics.ComputationalMode.Hankel;

mt = Numerics.ModalTruncation(H, contour, rd);
mt.compute();

H_in  = mt.getRegionTransferFunction();     % approximates H_D    (unstable part)
H_out = mt.getResidualTransferFunction();   % H - H_in            (stable complement)
ew    = mt.getRegionEigenvalues();          % eigenvalues recovered in D

%% Report recovered eigenvalues and decomposition accuracy
test_z   = 1i * logspace(-1, 1, 5);
pole_err = gmatch(ew_ust, ew);
dec_err  = arrayfun(@(z) abs(H(z) - (H_in(z) + H_out(z))), test_z);
fprintf('  Hankel: %d ew in D | pole err = %.3e | max|H-(H_D+H_Dc)| = %.3e\n', ...
        numel(ew), pole_err, max(dec_err));

%% Figure 1: contour + true vs. computed eigenvalues (zoomed near contour)
% The full Boeing spectrum sprawls (|Im| up to ~50, Re down to ~-540), so we
% lock the view tightly around the contour. Order matters: `axis equal` adjusts
% aspect but can expand limits to fit data; setting xlim/ylim *after* it locks
% the view to the zoomed window around the unstable cluster.
figure(1); clf; hold on;
t_arc = linspace(phi(1), phi(2), 300);
plot(real(gamma + rho*exp(1i*t_arc)), imag(gamma + rho*exp(1i*t_arc)), ...
     'k-', 'LineWidth', 1.5, 'DisplayName', 'contour');
plot(real(gamma + rho*[exp(1i*phi(2)), exp(1i*phi(1))]), ...
     imag(gamma + rho*[exp(1i*phi(2)), exp(1i*phi(1))]), ...
     'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
scatter(real(ew_stab), imag(ew_stab), 60, 'b', 'filled', ...
        'DisplayName', 'stable ($D^c$)');
scatter(real(ew_ust),  imag(ew_ust),  120, 'r', 'filled', ...
        'DisplayName', 'unstable ($D$)');
scatter(real(ew), imag(ew), 180, 'k', 'x', 'LineWidth', 2.5, ...
        'DisplayName', 'computed (Hankel)');
hold off; grid on;
axis equal;
% Pad ~ one contour radius beyond the disk on each side; the chord at Re=0
% means the left edge needs a sliver of negative-real margin to read clean.
xlim([-0.5*rho, 1.6*rho]); ylim([-1.5*rho, 1.5*rho]);
hold on; plot([0 0], ylim, 'k:', 'LineWidth', 0.5, 'HandleVisibility', 'off');
hold off;
xlabel('$\Re z$', 'Interpreter', 'latex'); ylabel('$\Im z$', 'Interpreter', 'latex');
legend('Location', 'northoutside', 'Orientation', 'horizontal', ...
       'Interpreter', 'latex');
ax = gca; ax.Toolbar = [];
fontsize(fsp, 'points');
if run_mode == "figs"
    exportgraphics(gcf, fullfile(fig_dir,'mt_boeing_poles.pdf'), ...
                   'ContentType','vector');
end

%% Figure 2: magnitudes of H, H_D, H_{D^c} across frequency
% Boeing low-frequency content sits well below 1 rad/s; sweep four decades.
wband  = logspace(-2, 2, 800);
H_mag  = zeros(size(wband)); HD_mag = H_mag; HDc_mag = H_mag;
for k = 1:numel(wband)
    zk = 1i*wband(k);
    H_mag(k)   = abs(H(zk));
    HD_mag(k)  = abs(H_in(zk));
    HDc_mag(k) = abs(H_out(zk));
end
figure(2); clf;
loglog(wband, H_mag,   'k-',  'LineWidth', 2.5, 'DisplayName', '$|H|$'); hold on;
loglog(wband, HD_mag,  'r--', 'LineWidth', 1.8, 'DisplayName', '$|\hat{H}_{D}|$');
loglog(wband, HDc_mag, 'b--', 'LineWidth', 1.8, 'DisplayName', '$|\hat{H}_{D^c}|$');
hold off; grid on;
xlabel('$\omega$', 'Interpreter', 'latex');
ylabel('magnitude $|H(\mathrm{i}\omega)|$', 'Interpreter', 'latex');
legend('Location', 'northoutside', 'Orientation', 'horizontal', ...
       'Interpreter', 'latex');
ax = gca; ax.Toolbar = [];
fontsize(fsp, 'points');
if run_mode == "figs"
    exportgraphics(gcf, fullfile(fig_dir,'mt_boeing_bode.pdf'), ...
                   'ContentType','vector');
end

if run_mode == "figs"
    fprintf('Figures written to %s\n', fig_dir);
end
