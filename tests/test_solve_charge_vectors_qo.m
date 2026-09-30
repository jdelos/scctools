function test_solve_charge_vectors_qo()
%TEST_SOLVE_CHARGE_VECTORS_QO Keep Qo type aligned with duty type.
if exist('OCTAVE_VERSION','builtin'), pkg load symbolic; end
addpath(fileparts(fileparts(mfilename('fullpath'))));
Q = {[-1 1 1; 1 -1 0], [-1 1 0; 1 -1 1]};
[al_numeric,m_numeric] = solve_charge_vectors(Q,1,0.5);
assert(isnumeric(al_numeric{1}) && isnumeric(m_numeric));
D = sym('D');
[al_symbolic,m_symbolic] = solve_charge_vectors(Q,1,D);
assert(isa(al_symbolic{1},'sym') && isa(m_symbolic,'sym'));
assert(isequal(size(al_symbolic{1}),size(al_numeric{1})));
assert(isequal(size(m_symbolic),size(m_numeric)));
end
