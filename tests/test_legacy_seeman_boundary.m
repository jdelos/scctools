function test_legacy_seeman_boundary()
%TEST_LEGACY_SEEMAN_BOUNDARY Compatibility entry point remains callable.
% Detailed schema and filtering coverage lives in test_generic_qfa_boundary.
addpath(fileparts(fileparts(mfilename('fullpath'))));
if exist('OCTAVE_VERSION','builtin'), pkg load symbolic; end
arch = dickson_arch(2);
result = hybrid_topology(arch, sym(1)/2, struct('dc_out', true));
assert(isstruct(result) && isfield(result, 'is'));
end
