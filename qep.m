% Figures land in code/figures/ next to this script. CIMTOOL is assumed to
% already be on the MATLAB path -- see code/README.md for setup.
thisfile = mfilename('fullpath');
if isempty(thisfile); thisfile = fullfile(pwd,'qep.m'); end
figDir = fullfile(fileparts(thisfile),'figures');
if ~exist(figDir,'dir'); mkdir(figDir); end

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

% CIMTOOL is an App Designer GUI -- skip when MATLAB was started with -batch
% (the figure backend isn't available, and the script otherwise hangs).
if usejava('desktop')
    cimtool = CIMTOOL(cim, vtTheme);
end

%% some text modes
cim.setComputationalMode(Numerics.ComputationalMode.Hankel);

% cim.SampleData.ell = 2; cim.SampleData.r = 2;
% cim.RealizationData.RealizationSize = Numerics.RealizationSize(4,2);
% cim.SampleData.Contour.N = 128;

cim.SampleData.ell = 2; cim.SampleData.r = 2;
cim.RealizationData.RealizationSize = Numerics.RealizationSize(4,8);
cim.SampleData.Contour.N = 64;

ew = cim.eigs(); ERA_err = gmatch([-1;1;1;1],ew)