# cim-ymmor-2026

Companion MATLAB code for the YMMOR 2026 talk
**Contour Integral Methods for the Masses** by Dan Folescu
(joint with Mark Embree and Serkan Gugercin, Virginia Tech).

The talk introduces [**CIMTOOL**](https://github.com/dan123222123/CIMTOOL), a
MATLAB package for solving linear / generalized / nonlinear eigenvalue problems
by contour integration, with realization via Hankel / ERA, SPLoewner (single-point Loewner),
and MPLoewner (multi-point Loewner).

## Assumptions

The drivers make two assumptions, and do nothing to discover or work
around either:

1. **CIMTOOL is already on the MATLAB path.** The scripts call
   `Numerics.ModalTruncation`, `Numerics.sploewner.*`, `Numerics.realize`,
   `Visual.Contour.Ellipse`, `Visual.OperatorData`, `Visual.SampleData`,
   `Visual.CIM`, and `Visual.CIMTOOL` directly. Put CIMTOOL on the path
   before running any driver, e.g. via your MATLAB `startup.m`:

   ```matlab
   addpath(genpath('/path/to/CIMTOOL/src'));
   ```

   or interactively / per-session:

   ```matlab
   addpath(genpath('/path/to/CIMTOOL/src'));
   ```

2. **Figures are written to `code/figures/`** (a sibling of the driver
   scripts). The folder is created on first run, is `.gitignore`'d, and
   is computed from `mfilename('fullpath')` so it does not depend on the
   current working directory. If you want figures elsewhere, edit
   `figDir = fullfile(...,'figures')` in the script's first cell.

## Quickstart

```matlab
% In MATLAB, with CIMTOOL on the path and this folder as cwd:
modal_truncation              % synthetic 2-RHP-pole / 4-LHP-pole split
modal_truncation_boeing767    % Boeing 767 SISO, real stable-unstable split
qep                           % Tisseur & Meerbergen QEP, CIMTOOL launcher
exact_data                    % exact-Markov-parameter \sigma sweep
quadrature_data               % same sweep with quadrature samples
```

End-to-end from a shell:

```bash
matlab -batch "addpath(genpath('/path/to/CIMTOOL/src')); \
               addpath('/path/to/code'); modal_truncation"
```

## Dependencies

| Component | Version tested | How to get it |
| --- | --- | --- |
| MATLAB | R2024b / R2025b | https://www.mathworks.com (any reasonably recent release should work; the scripts use `arguments`, `string`, and `exportgraphics`) |
| CIMTOOL | latest `main` | `git clone https://github.com/dan123222123/CIMTOOL` |
| Parallel Computing Toolbox | optional | only `exact_sploewner_sigma_choice.m` and `qep.m` use `parfor`; serial `for` works too |

No proprietary data is required: `data/boeing767_nnLTI.mat` ships with the repo.

## Layout

```
code/
├── README.md
├── .gitignore
├── data/
│   └── boeing767_nnLTI.mat              # Boeing 767 aeroelastic model (n=55)
├── figures/                             # driver outputs land here (created on first run, .gitignore'd)
├── modal_truncation.m                   # Demo 4a: synthetic stable/unstable split
├── modal_truncation_boeing767.m         # Demo 4b: Boeing 767 SISO, ELEV-QCG channel
├── qep.m                                # NEP example: gyroscopic QEP (Tisseur & Meerbergen 2001 §3.6)
├── exact_data.m                         # Where should sigma live? -- exact Markov parameters
├── quadrature_data.m                    # Same sweep with quadrature-sampled moments
├── gmatch.m                             # greedy 1-1 matching of computed -> reference eigenvalues
├── realize_inorder.m                    # eig(Ds, Db) sorted by descending |lambda|
├── redblue.m                            # diverging colormap (BSD, Auton 2009)
├── redblue_license.txt
└── html/                                # `matlab publish` output (not under version control by default)
```

## License

Code in this folder is released under the MIT License (see `LICENSE` once
added). `redblue.m` is third-party, BSD-licensed -- the original notice is
preserved in `redblue_license.txt`.

## Acknowledgements

This work is supported by **NSF DMS-2411141**. Joint with Mark Embree and
Serkan Gugercin (Virginia Tech). The CIMTOOL package and these drivers
build on prior contour-integral work by Beyn, Brennan, Embree, and others;
see references in the talk slides for details.
