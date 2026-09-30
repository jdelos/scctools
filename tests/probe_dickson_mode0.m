function report = probe_dickson_mode0()
%PROBE_DICKSON_MODE0 Probe public normal-mode Dickson construction.
% Failure reports preserve graph/tree/cut-set shapes; no production behavior changes.
if exist('OCTAVE_VERSION','builtin'), pkg load symbolic; end
addpath(fileparts(fileparts(mfilename('fullpath'))));
report = struct('schema_version','scctools.issue27.mode0.v2', ...
    'serializer_version','ordered-string-matrix.v1', ...
    'runtime',version,'package','symbolic','commit',git_commit(), ...
    'units','incidence entries dimensionless; duties fractions; A and m normalized', ...
    'assumptions','D symbolic; dc_out=true; half_point=false; normal mode means no Mode argument; exact shapes', ...
    'tolerances','exact symbolic equality; no numeric tolerance', 'error_cases',error_cases(), 'cases',[]);
for n_caps = [2 3]
    c = struct('n_caps',n_caps,'architecture',encode_arch(dickson_arch(n_caps)), ...
        'adapter',empty_result(), ...
        'direct_generic',empty_result(),'comparison',struct());
    c.adapter = run_adapter(n_caps);
    c.direct_generic = run_direct(n_caps);
    if c.adapter.success ~= c.direct_generic.success
        error('probe:Mismatch','Adapter/direct success mismatch for n_caps=%d',n_caps);
    end
    if c.adapter.success
        if ~isequal(c.adapter.m_ratios,c.direct_generic.m_ratios) || ...
                ~isequal(c.adapter.duties,c.direct_generic.duties) || ...
                ~isequal(c.adapter.A,c.direct_generic.A)
            error('probe:Mismatch','Adapter/direct serialized output mismatch for n_caps=%d',n_caps);
        end
        assert(isstruct(c.adapter.m_ratios) && isfield(c.adapter.m_ratios,'values'));
        c.comparison.m_ratios_equal = true;
        c.comparison.duties_equal = true;
        c.comparison.A_equal = true;
    else
        c.comparison.error_identifiers_equal = strcmp(c.adapter.error.identifier,c.direct_generic.error.identifier);
        c.comparison.error_messages_equal = strcmp(c.adapter.error.message,c.direct_generic.error.message);
        c.comparison.instrumentation = instrument_case(n_caps);
    end
    report.cases = [report.cases c];
end
fprintf('%s\n',jsonencode(report));
end

function r = empty_result()
r = struct('success',false,'m_ratios',[],'duties',[],'A',[],'graph',[],'tree_indices',[],'cutset',[], ...
    'branch_metadata',[],'phase_count',[],'ordered_symbols',[],'error',struct(), ...
    'provenance','generic_switched_capacitor_class(dickson_arch(n_caps),''Duty'',sym(''D''))');
end
function r = run_adapter(n)
r=empty_result();
r.provenance='dickson_hybrid_topology(n_caps,sym(''D''),struct(''dc_out'',true,''half_point'',false))';
try
    t=dickson_hybrid_topology(n,sym('D'),struct('dc_out',true,'half_point',false));
    assert(isstruct(t) && isfield(t,'g_top'),'adapter construction failed');
    r.success=true; r.m_ratios=encode_matrix(t.g_top.m_ratios); r.duties=encode_matrix(t.g_top.duty);
    r=characterize(t.g_top,r);
catch e
    r.error=err_struct(e);
end
end
function r = run_direct(n)
r=empty_result();
r.provenance='generic_switched_capacitor_class(dickson_arch(n_caps),''Duty'',sym(''D''))';
try
    t=generic_switched_capacitor_class(dickson_arch(n),'Duty',sym('D'));
    assert(isa(t,'generic_switched_capacitor_class'),'generic construction failed');
    r.success=true; r.m_ratios=encode_matrix(t.m_ratios); r.duties=encode_matrix(t.duty);
    r=characterize(t,r);
catch e
    r.error=err_struct(e);
end
end
function r = characterize(t,r)
r.phase_count=t.n_phases;
r.A=cell(1,t.n_phases); r.graph=cell(1,t.n_phases); r.tree_indices=cell(1,t.n_phases); r.cutset=cell(1,t.n_phases); r.branch_metadata=cell(1,t.n_phases);
for p=1:t.n_phases
    ph=t.phase{p}; g=ph.get_on_no_sw(); tree=tree_ph_scc(g,t.n_caps,0); q=fun_cutset(g,tree);
    r.A{p}=encode_matrix(ph.get_a_vector()); r.graph{p}=encode_matrix(g); r.tree_indices{p}=encode_matrix(tree); r.cutset{p}=encode_matrix(q);
    r.branch_metadata{p}=struct('inc_on_conv_sw',encode_matrix(ph.inc_on_conv_sw),'inc_on_conv',encode_matrix(ph.inc_on_conv),'sw_idxs',encode_matrix(ph.sw_idxs),'n_caps',ph.n_caps,'n_loads',ph.n_loads,'n_on_sw',ph.n_on_sw,'n_off_sw',ph.n_off_sw);
end
r.ordered_symbols=ordered_symbols(t);
end
function s = ordered_symbols(t)
s={};
for p=1:t.n_phases
    vals={t.phase{p}.inc_on_conv_sw,t.phase{p}.inc_on_conv};
    for k=1:numel(vals)
        v=symvar(vals{k});
        for j=1:numel(v), name=char(v(j)); if ~any(strcmp(s,name)), s{end+1}=name; end, end
    end
end
end
function x = encode_arch(a)
x=struct('Acaps',encode_matrix(a.Acaps),'Asw',encode_matrix(a.Asw), ...
    'Asw_act',encode_matrix(a.Asw_act));
end
function x = encode_matrix(a)
% Exact strings, flat row-major values, explicit shape/order metadata.
x=struct('rows',size(a,1),'cols',size(a,2),'order','row-major','values',{{}});
x.values=cell(1,numel(a));
k=0;
for i=1:size(a,1)
    for j=1:size(a,2), k=k+1; x.values{k}=encode_scalar(a(i,j)); end
end
end
function x = encode_scalar(v)
if isa(v,'sym'), x=char(v);
elseif isnumeric(v) || islogical(v), x=sprintf('%.17g',double(v));
else, error('probe:UnsupportedValue','Unsupported fixture value class %s',class(v)); end
end
function x = err_struct(e)
x=struct('identifier',e.identifier,'message',e.message,'stack',[]);
for k=1:numel(e.stack)
    x.stack=[x.stack struct('file',e.stack(k).file,'name',e.stack(k).name,'line',e.stack(k).line)];
end
end
function x = instrument_case(n)
x=struct('phase_graph_shapes',[],'tree_shapes',[],'cutset_shapes',[],'error',struct());
try
    arch=dickson_arch(n); nodes=size(arch.Acaps,1); loads=eye(nodes); loads(:,1)=[];
    supply=zeros(nodes,1); supply(1)=1; d=sym('D'); x.phase_graph_shapes=[];
    for p=1:2
        on=arch.Asw(:,logical(arch.Asw_act(p,:))); off=arch.Asw(:,~logical(arch.Asw_act(p,:)));
        ph=SCC_Phase(on,arch.Acaps,off,loads*(d*(1-0)+0),size(arch.Acaps,2),find(arch.Asw_act(p,:)),supply);
        g=ph.get_on_no_sw(); t=tree_ph_scc(g,size(arch.Acaps,2),0); q=fun_cutset(g,t);
        x.phase_graph_shapes=[x.phase_graph_shapes struct('graph',size(g),'tree',size(t),'cutset',size(q))];
    end
catch e
    x.error=err_struct(e);
end
end
function x=error_cases()
% Public adapter has no documented invalid-architecture or phase-count API.
x=struct('invalid_architecture',struct('status','not-representable','metadata','No public adapter entry point accepts an arbitrary architecture.'), ...
    'singular_graph',struct('status','not-representable','metadata','No public adapter entry point exposes singular graph construction.'), ...
    'more_than_two_phases',struct('status','not-representable','metadata','dickson_hybrid_topology public adapter constructs fixed two-phase architecture; no semantic-changing probe.'));
end
function c=git_commit()
[~,c]=system('git rev-parse HEAD'); c=strtrim(c);
end
