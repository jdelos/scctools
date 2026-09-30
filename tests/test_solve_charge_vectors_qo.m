function test_solve_charge_vectors_qo()
%TEST_SOLVE_CHARGE_VECTORS_QO Keep Qo type aligned with duty type.
if exist('OCTAVE_VERSION','builtin'), pkg load symbolic; end
addpath(fileparts(fileparts(mfilename('fullpath'))));
Q_numeric = {[-1 1 2; 1 0 3], [-1 0 1/2; 0 1 2]};
[al_numeric,m_numeric] = solve_charge_vectors(Q_numeric,1,0.5);
assert(isnumeric(al_numeric{1}) && isnumeric(m_numeric));
assert(isequal(size(al_numeric{1}),[2 1]));
assert(isequal(size(al_numeric{2}),[2 1]));
assert(isequal(size(m_numeric),[1 1]));
assert(isequal(al_numeric{1},[3;-5]));
assert(isequal(al_numeric{2},[-1/2;3]));
assert(isequal(m_numeric,2.5));

D = sym('D');
Q_symbolic = {Q_numeric{1}, [-1 0 D; 0 1 4*D]};
[al_symbolic,m_symbolic] = solve_charge_vectors(Q_symbolic,1,D);
assert(isa(al_symbolic{1},'sym') && isa(m_symbolic,'sym'));
assert(isequal(size(al_symbolic{1}),[2 1]));
assert(isequal(size(al_symbolic{2}),[2 1]));
assert(isequal(size(m_symbolic),[1 1]));
assert(isequal(al_symbolic{1},sym([3;-5])));
assert(isequal(al_symbolic{2},[-D;5-4*D]));
assert(isequal(m_symbolic,3-D));
end
