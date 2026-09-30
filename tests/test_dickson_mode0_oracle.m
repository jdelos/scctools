function test_dickson_mode0_oracle()
%TEST_DICKSON_MODE0_ORACLE Compare generated report with MATLAB R2021a fixture.
addpath(fileparts(fileparts(mfilename('fullpath'))));
fixture_path = fullfile(fileparts(mfilename('fullpath')),'fixtures','dickson_mode0_matlab_r2021a.json');
expected = jsondecode(fileread(fixture_path));
actual = probe_dickson_mode0();
assert(strcmp(actual.schema_version,'scctools.issue27.mode0.v2'));
assert(strcmp(actual.serializer_version,'ordered-string-matrix.v1'));
assert(~isempty(actual.commit) && ~isempty(regexp(actual.commit,'^[0-9a-f]{40}$','once')));
assert(~isempty(actual.runtime));
assert(~isempty(regexp(actual.runtime,'R2021a','once')) || ~isempty(regexp(actual.runtime,'^9\\.10\\.', 'once')));
assert(numel(actual.cases) == 2);
for k = 1:numel(actual.cases)
    c = actual.cases(k);
    assert(c.adapter.success && c.direct_generic.success);
    assert(strcmp(c.adapter.provenance,'dickson_hybrid_topology(n_caps,sym(''D''),struct(''dc_out'',true,''half_point'',false))'));
    assert(strcmp(c.direct_generic.provenance,'generic_switched_capacitor_class(dickson_arch(n_caps),''Duty'',sym(''D''))'));
    assert_serialized_architecture(c.architecture);
    assert_serialized_result(c.adapter);
    assert_serialized_result(c.direct_generic);
    assert_supplemental(c.adapter);
    assert_supplemental(c.direct_generic);
end
assert(isfield(actual,'error_cases'));
assert(all(isfield(actual.error_cases,{'invalid_architecture','singular_graph','more_than_two_phases'})));
actual = rmfield(actual,{'runtime','commit'});
% JSON round-trip normalizes MATLAB cell/struct container differences without
% changing flat row-major values or provenance strings.
actual = jsondecode(jsonencode(actual));
expected = jsondecode(jsonencode(expected));
expected = rmfield(expected,{'runtime','commit'});
assert(isequal(actual,expected),'Generated report differs from checked-in oracle');
end
function assert_serialized_architecture(a)
assert_serialized_matrix(a.Acaps); assert_serialized_matrix(a.Asw); assert_serialized_matrix(a.Asw_act);
end
function assert_serialized_result(r)
assert(r.success); assert_serialized_matrix(r.m_ratios); assert_serialized_matrix(r.duties);
for k = 1:numel(r.A), assert_serialized_matrix(r.A{k}); end
end
function assert_supplemental(r)
assert(isfield(r,'graph') && isfield(r,'tree_indices') && isfield(r,'cutset'));
assert(isfield(r,'branch_metadata') && isfield(r,'phase_count') && isfield(r,'ordered_symbols'));
assert(r.phase_count == numel(r.graph) && r.phase_count == numel(r.tree_indices));
assert(r.phase_count == numel(r.cutset) && r.phase_count == numel(r.branch_metadata));
for k = 1:r.phase_count
    assert_serialized_matrix(r.graph{k}); assert_serialized_matrix(r.tree_indices{k});
    assert_serialized_matrix(r.cutset{k});
    b = r.branch_metadata{k};
    assert(all(isfield(b,{'inc_on_conv_sw','inc_on_conv','sw_idxs','n_caps','n_loads','n_on_sw','n_off_sw'})));
    assert_serialized_matrix(b.inc_on_conv_sw); assert_serialized_matrix(b.inc_on_conv);
    assert_serialized_matrix(b.sw_idxs);
end
assert(iscell(r.ordered_symbols));
end
function assert_serialized_matrix(x)
assert(isfield(x,'rows') && isfield(x,'cols') && isfield(x,'order') && isfield(x,'values'));
assert(strcmp(x.order,'row-major')); assert(numel(x.values) == x.rows*x.cols);
end
