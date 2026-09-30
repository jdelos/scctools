function test_dickson_hybrid_boundary()
%TEST_DICKSON_HYBRID_BOUNDARY Check public defaults and supported options.
addpath(fileparts(fileparts(mfilename('fullpath'))));
if exist('OCTAVE_VERSION','builtin'), pkg load symbolic; end

% One-argument and two-argument calls retain documented default duty.
t_default = dickson_hybrid_topology(2);
t_two_arg = dickson_hybrid_topology(2, []);
assert(isequal(t_default.duty, 0.5));
assert(isequal(t_two_arg.duty, 0.5));

% Supported options remain usable.
t_no_dc = dickson_hybrid_topology(2, sym(1)/2, struct('dc_out', false, 'half_point', false));
assert(isstruct(t_no_dc) && isfield(t_no_dc, 'ratio'));

% Representable invalid architecture argument fails with stable identifier.
try
    dickson_hybrid_topology(0);
    error('test_dickson_hybrid_boundary:ExpectedFailure', 'n_caps=0 was accepted');
catch e
    assert(strcmp(e.identifier, 'dickson_hybrid_topology:InvalidNCaps'));
end

% Unsupported half-point mode fails before undefined matrices are referenced.
try
    dickson_hybrid_topology(2, sym(1)/2, struct('half_point', true));
    error('test_dickson_hybrid_boundary:ExpectedFailure', 'half_point=true was accepted');
catch e
    assert(strcmp(e.identifier, 'dickson_hybrid_topology:HalfPointUnsupported'));
end
end
