function test_fun_loop_zero_rows()
%TEST_FUN_LOOP_ZERO_ROWS Symbolic row deletion must use Octave-compatible indices.
if exist('OCTAVE_VERSION','builtin'), pkg load symbolic; end
addpath(fileparts(fileparts(mfilename('fullpath'))));
A = sym([0 0 0; 1 0 0]);
[Bf,T] = fun_loop(A,1);
assert(isequal(T,1));
assert(isequal(Bf,sym([0 1 0; 0 0 1])));
end
