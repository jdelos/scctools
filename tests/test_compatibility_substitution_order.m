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
    canonical_c = sym('C', [1 numel(t.var_ssl)]);
    canonical_ron = sym('Ron', [1 numel(t.var_fsl)]);
    canonical_cesr = sym('Resr', [1 numel(t.var_fesr)]);
    c_values = 2:2 + numel(canonical_c) - 1;
    ron_values = 11:11 + numel(canonical_ron) - 1;
    cesr_values = 21:21 + numel(canonical_cesr) - 1;
    assert(isequal(t.var_ssl, canonical_c));
    assert(isequal(t.var_fsl, canonical_ron));
    assert(isequal(t.var_fesr, canonical_cesr));
    assert(isequal(t.r_vars, canonical_c));
    assert(all(isAlways(t.eval_ssl(c_values) == subs(t.f_ssl, canonical_c, c_values)), 'all'));
    assert(all(isAlways(t.eval_fsl(ron_values) == subs(t.f_fsl, canonical_ron, ron_values)), 'all'));
    assert(all(isAlways(t.eval_fesr(cesr_values) == subs(t.f_esr, canonical_cesr, cesr_values)), 'all'));
    assert(all(isAlways(t.eval_r(c_values) == subs(t.r, canonical_c, c_values)), 'all'));
    assert(all(isAlways(t.eval_q_dc(c_values) == subs(t.q_dc, canonical_c, c_values)), 'all'));
end
end
