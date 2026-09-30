function test_legacy_seeman_boundary()
%TEST_LEGACY_SEEMAN_BOUNDARY Exercise actual Seeman producer and consumer.
addpath(fileparts(fileparts(mfilename('fullpath'))));
if exist('OCTAVE_VERSION','builtin'), pkg load symbolic; end
seeman = generate_topology('Dickson', 2);
assert(strcmp(seeman.schema, 'scctools.matlab.seeman.v1'));
assert(strcmp(seeman.topName, 'Dickson'));
techlib;
result = implement_topology(seeman, 1, flatMOS, flatcap);
assert(isstruct(result) && isfield(result, 'topology'));
try
    implement_topology(hybrid_topology(dickson_arch(2)), 1, flatMOS, flatcap);
    error('test_legacy_seeman_boundary:MissingSchemaGuard');
catch err
    assert(strcmp(err.identifier, 'implement_topology:InvalidSchema'));
end
end
