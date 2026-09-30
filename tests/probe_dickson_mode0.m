function report = probe_dickson_mode0()
%PROBE_DICKSON_MODE0 Probe public normal-mode Dickson construction.
% Failure reports preserve graph/tree/cut-set shapes; no production behavior changes.
if exist('OCTAVE_VERSION','builtin'), pkg load symbolic; end
addpath(fileparts(fileparts(mfilename('fullpath'))));
report = struct('schema_version','scctools.issue27.mode0.v1', ...
    'runtime',version,'package','symbolic','commit',git_commit(), ...
    'units','incidence entries dimensionless; duties fractions; A and m normalized', ...
    'assumptions','D symbolic; dc_out=true; half_point=false; normal mode means no Mode argument; exact shapes', ...
    'tolerances','exact symbolic equality; no numeric tolerance', 'cases',[]);
for n_caps = [2 3]
    c = struct('n_caps',n_caps,'adapter',empty_result(), ...
        'direct_generic',empty_result(),'comparison',struct());
    c.adapter = run_adapter(n_caps);
    c.direct_generic = run_direct(n_caps);
    if c.adapter.success ~= c.direct_generic.success
        error('probe:Mismatch','Adapter/direct success mismatch for n_caps=%d',n_caps);
    end
    if c.adapter.success
        if ~isequal(c.adapter.m_ratios,c.direct_generic.m_ratios)
            error('probe:Mismatch','Adapter/direct m_ratios mismatch for n_caps=%d',n_caps);
        end
        c.comparison.m_ratios_equal = true;
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
r = struct('success',false,'m_ratios',[],'duties',[],'A',[],'error',struct(), ...
    'provenance','dickson_hybrid_topology(n_caps,sym(''D''),struct(''dc_out'',true,''half_point'',false))');
end
function r = run_adapter(n)
r=empty_result();
try
    t=dickson_hybrid_topology(n,sym('D'),struct('dc_out',true,'half_point',false));
    assert(isstruct(t) && isfield(t,'g_top'),'adapter construction failed');
    r.success=true; r.m_ratios=t.g_top.m_ratios; r.duties=t.g_top.duty;
    r.A=cell(1,numel(t.g_top.phase));
    for p=1:numel(t.g_top.phase), r.A{p}=t.g_top.phase{p}.get_a_vector(); end
catch e
    r.error=err_struct(e);
end
end
function r = run_direct(n)
r=empty_result();
try
    t=generic_switched_capacitor_class(dickson_arch(n),'Duty',sym('D'));
    assert(isa(t,'generic_switched_capacitor_class'),'generic construction failed');
    r.success=true; r.m_ratios=t.m_ratios; r.duties=t.duty;
    r.A=cell(1,numel(t.phase));
    for p=1:numel(t.phase), r.A{p}=t.phase{p}.get_a_vector(); end
catch e
    r.error=err_struct(e);
end
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
function c=git_commit()
[~,c]=system('git rev-parse HEAD'); c=strtrim(c);
end
