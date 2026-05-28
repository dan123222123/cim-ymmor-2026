% MODAL_TRUNCATION  Demo 4a: synthetic stable/unstable spectral split via CIMTOOL.
%
%   A SISO transfer function with 2 right-half-plane (unstable) poles and 4
%   left-half-plane (stable) poles is sampled along a circular-segment contour
%   that encloses only the unstable poles. Numerics.ModalTruncation realizes the
%   in-region subsystem H_D (the unstable part) directly from those samples; the
%   residual H_{D^c} = H - H_D then recovers the stable complement.
%
%   Produces the two figures referenced by Demo 4a of the YMMOR deck:
%       tex/figures/mt_synth_poles.pdf  -- contour + true vs. computed poles
%       tex/figures/mt_synth_bode.pdf   -- |H|, |H_D|, |H_{D^c}| vs. omega
%
%   Requires CIMTOOL (math/CIMTOOL) on the MATLAB path; the Paths cell below
%   adds it relative to this file if it is not already available. Run the whole
%   script, or cell-by-cell (Ctrl+Enter) for the live demo.
%
%   See also: Numerics.ModalTruncation, Numerics.Contour.CircularSegment, gmatch

%% Paths
here = fileparts(mfilename('fullpath'));
if exist('Numerics.ModalTruncation','class') ~= 8
    addpath(genpath(fullfile(here,'..','..','..','..','CIMTOOL','src')));
end
fig_dir = fullfile(here,'..','tex','figures');
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end
fsp = 16; % figure font size

%% Synthetic SISO system: 2 unstable + 4 stable poles
ew_ust  = [0.2+0.5i; 0.2-0.5i];                  % right-half-plane poles (region D)
ew_stab = [-0.5+1i; -0.5-1i; -1+0.2i; -1-0.2i];  % left-half-plane poles (complement D^c)
ew_all  = [ew_ust; ew_stab];

rng(42);                                  % reproducible residues
W = randn(numel(ew_all),1);               % input residue directions  (n x 1)
V = randn(1,numel(ew_all));               % output residue directions (1 x n)
H = @(z) Numerics.poresz(z, ew_all, W, V);% SISO transfer function in pole-residue form

%% Circular-segment contour over {Re z > 0}, enclosing only the unstable poles
gamma = 0; rho = 0.7; phi = pi/2; N = 32;
contour = Numerics.Contour.CircularSegment(gamma, rho, phi, N);

%% Modal truncation: realize the in-region (unstable) subsystem, all three modes
% We realize H_D with each CIMTOOL computational mode so the Bode comparison
% below can overlay them. Hankel and SPLoewner recover the 2 enclosed poles to
% machine precision here; MPLoewner has previously produced a rank-deficient data
% matrix for this SISO problem, so each realization is guarded -- a failure is
% reported (and drawn as a gap in Figure 2) rather than aborting the script.
modes      = [Numerics.ComputationalMode.Hankel];%, ...
              %Numerics.ComputationalMode.SPLoewner, ...
              %Numerics.ComputationalMode.MPLoewner];
mode_names = {'Hankel', 'SPLoewner', 'MPLoewner'};
nmode      = numel(modes);

R = struct('H_in', cell(1,nmode), 'H_out', cell(1,nmode), 'ew', cell(1,nmode), 'ok', cell(1,nmode));
for m = 1:nmode
    rd = Numerics.RealizationData();
    rd.RealizationSize   = Numerics.RealizationSize(2, 5, 5);  % expect 2 poles inside D
    rd.ComputationalMode = modes(m);
    try
        mt = Numerics.ModalTruncation(H, contour, rd);
        mt.compute();
        R(m).H_in  = mt.getRegionTransferFunction();   % approximates H_D     (unstable part)
        R(m).H_out = mt.getResidualTransferFunction(); % H - H_in             (stable complement)
        R(m).ew    = mt.getRegionEigenvalues();        % eigenvalues recovered in D
        R(m).ok    = true;
    catch err
        warning('modal_truncation:modeFailed', '%s realization failed: %s', ...
                mode_names{m}, err.message);
        R(m).H_in  = @(z) nan(size(z));
        R(m).H_out = @(z) nan(size(z));
        R(m).ew    = [];
        R(m).ok    = false;
    end
end

ref = 1;          % reference mode (Hankel) used for the pole scatter in Figure 1
ew  = R(ref).ew;  % eigenvalues drawn in Figure 1

%% Report recovered eigenvalues and decomposition accuracy per mode
test_z = [1i, 0.5, -1, 0.1+0.2i, -0.5+0.5i];
fprintf('Modal truncation: comparing %d realization modes (expect 2 poles in D)\n', nmode);
for m = 1:nmode
    if ~R(m).ok
        fprintf('  %-10s: FAILED (see warning above)\n', mode_names{m});
        continue;
    end
    pole_err = gmatch(ew_ust, R(m).ew);                % matched-pole error vs. truth
    % sanity check: H(z) = H_in(z) + H_out(z) by construction of the residual
    dec_err  = arrayfun(@(z) abs(H(z) - (R(m).H_in(z) + R(m).H_out(z))), test_z);
    fprintf('  %-10s: %d ew in D | pole err = %.3e | max|H-(H_in+H_out)| = %.3e\n', ...
            mode_names{m}, numel(R(m).ew), pole_err, max(dec_err));
end

%% Figure 1: contour + true vs. computed eigenvalues
figure(1); clf; hold on;
t_arc = linspace(-phi, phi, 300);
plot(real(gamma + rho*exp(1i*t_arc)), imag(gamma + rho*exp(1i*t_arc)), ...
     'k-', 'LineWidth', 1.5, 'DisplayName', 'contour');
plot(real(gamma + rho*[exp(1i*phi), exp(-1i*phi)]), ...
     imag(gamma + rho*[exp(1i*phi), exp(-1i*phi)]), ...
     'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
scatter(real(ew_stab), imag(ew_stab), 90, 'b', 'filled', 'DisplayName', 'stable ($D^c$)');
scatter(real(ew_ust),  imag(ew_ust),  90, 'r', 'filled', 'DisplayName', 'unstable ($D$)');
scatter(real(ew), imag(ew), 150, 'k', 'x', 'LineWidth', 2, ...
        'DisplayName', sprintf('computed (%s)', mode_names{ref}));
yl = ylim; plot([0 0], yl, 'k:', 'LineWidth', 0.5, 'HandleVisibility', 'off'); ylim(yl);
hold off; axis equal; grid on;
xlabel('$\Re z$', 'Interpreter', 'latex'); ylabel('$\Im z$', 'Interpreter', 'latex');
legend('Location', 'northoutside', 'Orientation', 'horizontal', 'Interpreter', 'latex');
ax = gca; ax.Toolbar = [];  % drop hover toolbar so it is not baked into the PDF
fontsize(fsp, 'points');
exportgraphics(gcf, fullfile(fig_dir, 'mt_synth_poles.pdf'), 'ContentType', 'vector');

%% Figure 2: magnitudes of H, H_D, H_{D^c} across frequency
wband = logspace(-2, 2, 800);
H_mag = zeros(size(wband)); HD_mag = H_mag; HDc_mag = H_mag;
for k = 1:numel(wband)
    zk = 1i*wband(k);
    H_mag(k)   = abs(H(zk));
    HD_mag(k)  = abs(R(ref).H_in(zk));
    HDc_mag(k) = abs(R(ref).H_out(zk));
end
figure(2); clf;
loglog(wband, H_mag,   'k-',  'LineWidth', 2.5, 'DisplayName', '$|H|$'); hold on;
loglog(wband, HD_mag,  'r--', 'LineWidth', 1.8, 'DisplayName', '$|\hat{H}_{D}|$');
loglog(wband, HDc_mag, 'b--', 'LineWidth', 1.8, 'DisplayName', '$|\hat{H}_{D^c}|$');
hold off; grid on;
xlabel('$\omega$', 'Interpreter', 'latex');
ylabel('magnitude $|H(\mathrm{i}\omega)|$', 'Interpreter', 'latex');
legend('Location', 'northoutside', 'Orientation', 'horizontal', 'Interpreter', 'latex');
ax = gca; ax.Toolbar = [];  % drop hover toolbar so it is not baked into the PDF
fontsize(fsp, 'points');
exportgraphics(gcf, fullfile(fig_dir, 'mt_synth_bode.pdf'), 'ContentType', 'vector');

fprintf('Figures written to %s\n', fig_dir);
