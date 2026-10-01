function test_qfa_bundle()
% Live native boundary vs MATLAB algebra; no expression-text comparisons.
request_file=[tempname '.json']; result_file=[tempname '.json'];
cleanup=onCleanup(@() remove_files(request_file,result_file)); %#ok<NASGU>
for nc=2:3
    req=struct('version',1,'architecture','qfa-graph','stages',2,'phases',2,'capacitors',nc,'duty',0.5);
    fprintf('QFA reference: %d capacitors\n',nc);
    ref=generic_switched_capacitor_class(dickson_arch(nc),'Duty',sym('D'));
    values=struct('duty',[0.25 0.75],'capacitances',1:nc,'switch_resistances',(1:ref.n_switches)/10,'capacitor_esr',(1:nc)/100,'frequency',{{7}});
    req.substitutions=values;
    fid=fopen(request_file,'w'); fprintf(fid,'%s\n',jsonencode(req)); fclose(fid);
    [status,msg]=system(sprintf('env -u LD_LIBRARY_PATH /tmp/scctools-qfa-json < "%s" > "%s"',request_file,result_file)); assert(status==0,msg);
    got=jsondecode(fileread(result_file)); assert(strcmp(got.type,'result'),jsonencode(got));
    fsw=sym('fsw');
    same(got.m,ref.m_ratios);
    for p=1:2
        same(phase(got.A,p),ref.phase{p}.get_a_vector());
        same(phase(got.B,p),ref.phase{p}.get_b_vector());
        same(phase(got.r,p),ref.phase{p}.get_r_vector());
        same(phase(got.G,p),ref.phase{p}.get_r_vector());
        same(phase(got.Ar,p),ref.phase{p}.get_ar_vector());
    end
    fprintf('QFA charge algebra passed: %d capacitors\n',nc);
    esr=subs(ref.k_fsl,ref.ron_switches,zeros(1,ref.n_switches));
    combined=sqrt((ref.k_ssl/fsw).^2+ref.k_fsl.^2);
    if nc==2
        same(got.ZSSL,ref.k_ssl/fsw); same(got.ZFSL,ref.k_fsl);
        same(got.ZESR,esr); same_squared(got.ZSCC,combined);
    else
        % ponytail: R2021a crashes expanding full three-cap loss identities.
        % Use exact unequal-component substitutions; restore full identities
        % when symbolic engine can finish them within bounded resources.
        ordered=[sym('D') ref.caps ref.ron_switches ref.esr_caps fsw];
        refs={ref.k_ssl/fsw,ref.k_fsl,esr,combined};
        names={'ZSSL','ZFSL','ZESR','ZSCC'};
        for sample=1:4
            replacements=[sym(sample)/5 sym((1:nc)+sample)/7 ...
                sym((1:ref.n_switches)+sample)/11 sym((1:nc)+sample)/101 sym(13+sample)];
            for f=1:4
                actual=subs(parse(got.(names{f})),ordered,replacements);
                expected=subs(refs{f},ordered,replacements);
                assert(isequal(size(actual),size(expected)));
                if f<4
                    assert_rational_zero(actual-expected);
                else
                    delta=vpa(actual-expected,50);
                    assert(all(isAlways(abs(delta(:))<sym(10)^(-40))),'RSS parity mismatch');
                end
            end
        end
    end
    fprintf('QFA loss algebra passed: %d capacitors\n',nc);
    same(got.stress.capacitor_voltage,ref.v_caps_norm);
    sw=sym(zeros(2,ref.n_switches)); current=sym(zeros(ref.n_switches,ref.n_outs));
    for p=1:2
        off=setdiff(1:ref.n_switches,ref.phase{p}.sw_idxs);
        sw(p,off)=ref.v_sw_norm(off);
        ar=ref.phase{p}.get_ar_vector(); current(ref.phase{p}.sw_idxs,:)=ar(nc+1:end,:);
    end
    same(got.stress.switch_voltage_by_phase,sw); same(got.stress.switch_current_by_output,current);
    assert(isequal(got.ordering.outputs(:),(2:ref.n_nodes)'));
    assert(isequal(got.symbols.capacitances(:),arrayfun(@char,ref.caps(:),'UniformOutput',false)));
    syms_order=[sym('D') ref.caps ref.ron_switches ref.esr_caps fsw];
    vals=[values.duty(1) values.capacitances values.switch_resistances values.capacitor_esr values.frequency{1}];
    references={ref.k_ssl/fsw,ref.k_fsl,esr,combined};
    fields={'ZSSL','ZFSL','ZESR','ZSCC'};
    for f=1:numel(fields)
        name=fields{f}; expected=double(subs(references{f},syms_order,vals));
        actual=double(parse(got.evaluated.(name)));
        assert(all(abs(actual(:)-expected(:))<1e-10));
    end
end
end
function x=phase(v,p)
if iscell(v) && ~ischar(v{1}) && ~isstring(v{1}), x=v{p}; else, x=squeeze(v(p,:,:)); end
end
function x=parse(v)
if iscell(v) && ~ischar(v{1}) && ~isstring(v{1})
    rows=cellfun(@(r) reshape(parse(r),1,[]),v,'UniformOutput',false); x=vertcat(rows{:});
elseif iscell(v)
    x=sym(zeros(size(v))); for k=1:numel(v), x(k)=str2sym(strrep(v{k},'**','^')); end
else, x=str2sym(strrep(v,'**','^')); end
end
function same(got,expected)
actual=parse(got); assert(isequal(size(actual),size(expected)),sprintf('QFA shape mismatch: %s vs %s',mat2str(size(actual)),mat2str(size(expected))));
assert_rational_zero(actual-expected);
end
function same_squared(got,expected)
actual=parse(got); assert(isequal(size(actual),size(expected)));
assert_rational_zero(actual.^2-expected.^2);
end
function assert_rational_zero(delta)
% Polynomial numerator identity avoids costly general rational factoring.
for k=1:numel(delta)
    [numerator,~]=numden(delta(k));
    assert(isequal(expand(numerator),sym(0)),'QFA algebra mismatch');
end
end
function remove_files(varargin)
for i=1:nargin, if isfile(varargin{i}), delete(varargin{i}); end, end
end
