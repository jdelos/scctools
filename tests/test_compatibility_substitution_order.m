function test_compatibility_substitution_order()
%TEST_COMPATIBILITY_SUBSTITUTION_ORDER Check positional variable order.
addpath(fileparts(fileparts(mfilename('fullpath'))));
if exist('OCTAVE_VERSION','builtin'), pkg load symbolic; end
adapters = { ...
    hybrid_topology(dickson_arch(2), sym(1)/2, struct('dc_out', true)), ...
    dickson_hybrid_topology(2, sym(1)/2, struct('dc_out', true, 'half_point', false)), ...
    ladder_hybrid_topology(2, sym(1)/2, struct('dc_out', true, 'half_point', false)), ...
    scc11_topology(sym(1)/2)};
for k = 1:numel(adapters)
    t = adapters{k};
    c = sym(2:2 + numel(t.var_ssl) - 1);
    ron = sym(11:11 + numel(t.var_fsl) - 1);
    cesr = sym(21:21 + numel(t.var_fesr) - 1);
    assert(all(isAlways(t.eval_ssl(c) == subs(t.f_ssl, t.var_ssl, c)), 'all'));
    assert(all(isAlways(t.eval_fsl(ron) == subs(t.f_fsl, t.var_fsl, ron)), 'all'));
    assert(all(isAlways(t.eval_fesr(cesr) == subs(t.f_esr, t.var_fesr, cesr)), 'all'));
    assert(all(isAlways(t.eval_r(c) == subs(t.r, t.r_vars, c)), 'all'));
    assert(all(isAlways(t.eval_q_dc(c) == subs(t.q_dc, t.r_vars, c)), 'all'));
    assert(numel(unique(c)) == numel(c));
    assert(numel(unique(ron)) == numel(ron));
    assert(numel(unique(cesr)) == numel(cesr));
end
end
