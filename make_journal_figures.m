% Makes figures for ReactiveShearInstability-Manuscript
addpath(genpath(pwd))
load_figure_options_journal
load_colormap
exDir=exist('Output','dir');
if ~exDir
    mkdir('Output',pwd)
end
% Figure 1 - Diagram (not created here)

% Figure 2 - Showing dependence of growth rate on angle of wavenumber
fig_wavenumber_angle

% Figure 3 - Showing dependence of growth rate on angle of wavenumber
% restricted to x-z plane
fig_wavenumber_angle_xz_plane

% Figure 4 - Showing dependence of growth rate on angle of wavenumber
% projected on kx-ky plane
fig_wavenumber_angle_kx_ky

% Figure 5 - Showing dependence of growth rate on compaction length
fig_compaction_length

% Figure 6 - MOR Diagram (not created here)

% Figure 7 - Showing background MOR state
fig_MOR_base_state

% Figure 8 - Showing evolution of wavenumber vector
fig_MOR_wavenumber

%Figure 9 - Showing upper bound on growth due to reaction and shear
fig_MOR_growth_upper

% Figure 10 - Showing effect of initial wavenumber orientation on growth
fig_MOR_growth_combined

% Figure 11 - Showing effect of a spectrum of initial wavenumber
% orientation
fig_MOR_growth_spectrum

% Figure 12 - Showing evolution of preferred orientation
fig_MOR_growth_favoured

% Figure 13 - Showing sensitivity (parameter choice)
fig_MOR_growth_favoured_sensitivity

% Figure 14 - Showing sensitivity (lithospheric dip)
fig_MOR_growth_combined_shallow_dip

% Figure 15 - Showing summary for discussion section
fig_MOR_summary

% Figure A1 - Showing alternative wavevector maximization procedure
fig_MOR_growth_local_max