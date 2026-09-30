function test_generic_qfa_boundary()
%TEST_GENERIC_QFA_BOUNDARY Keep generic QFA adapter schema explicit.
addpath(fileparts(fileparts(mfilename('fullpath'))));
if exist('OCTAVE_VERSION','builtin'), pkg load symbolic; end
fixture = jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')), 'fixtures', ...
    'hybrid_topology_compatibility.json')));
arch2 = dickson_arch(2);
arch3 = dickson_arch(3);

% MATLAB defaults remain callable.
defaults = hybrid_topology(arch2);
explicit = hybrid_topology(arch2, 0.5, struct('dc_out', true));
assert(strcmp(defaults.schema, fixture.legacy_schema));
assert(isequal(defaults.duty, 0.5));
assert(isequal(defaults.N_outs, explicit.N_outs));
assert(isequal(defaults.ratio, explicit.ratio));
assert(all(isfield(defaults, {'ratio','vc','vr','is','Y_ssl','Y_fsl','f_ssl', ...
    'f_fsl','f_esr','duty','g','dc_outputs','r','r_vars','q_dc','N_outs', ...
    'N_sw','N_caps','ph','g_top','schema'})));

% dc_out filtering stays explicit for both supported capacitor counts.
no_dc2 = hybrid_topology(arch2, 0.5, struct('dc_out', false));
no_dc3 = hybrid_topology(arch3, 0.5, struct('dc_out', false));
assert(no_dc2.N_outs == fixture.dc_out_false_n_outs_n2);
assert(no_dc3.N_outs == fixture.dc_out_false_n_outs_n3);
assert(no_dc2.N_outs < explicit.N_outs);
assert(no_dc3.N_outs < hybrid_topology(arch3, 0.5, struct('dc_out', true)).N_outs);

% Graph QFA and generic compatibility adapters cannot be mistaken for old schema.
modern = dickson_hybrid_topology(2, sym(1)/2, struct('dc_out', true, 'half_point', false));
assert(strcmp(modern.schema, fixture.graph_qfa_schema));
assert(~strcmp(defaults.schema, modern.schema));
end
