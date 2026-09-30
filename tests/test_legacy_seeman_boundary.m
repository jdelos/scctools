function test_legacy_seeman_boundary()
%TEST_LEGACY_SEEMAN_BOUNDARY Keep legacy Seeman topology callable separately.
addpath(fileparts(fileparts(mfilename('fullpath'))));
if exist('OCTAVE_VERSION','builtin'), pkg load symbolic; end

arch = dickson_arch(2);
legacy = hybrid_topology(arch, sym(1)/2, struct('dc_out', true));
assert(isstruct(legacy));
assert(all(isfield(legacy, {'ratio','vc','vr','is','Y_ssl','Y_fsl','f_ssl', ...
    'f_fsl','f_esr','duty','g','dc_outputs','r','r_vars','q_dc','N_outs', ...
    'N_sw','N_caps','ph','g_top'})));
assert(isequal(legacy.N_caps, 2));
assert(isequal(legacy.duty, sym(1)/2));

% Legacy and graph-QFA schemas stay distinct; adapter does not leak native fields.
modern = dickson_hybrid_topology(2, sym(1)/2, struct('dc_out', true, 'half_point', false));
assert(isstruct(modern) && isfield(modern, 'g_top'));
assert(~isfield(modern, 'phase_count'));
end
