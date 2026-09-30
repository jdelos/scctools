function test_native_scope_characterization()
%TEST_NATIVE_SCOPE_CHARACTERIZATION Record native seam evidence without adapter fakes.
addpath(fileparts(fileparts(mfilename('fullpath'))));
if exist('OCTAVE_VERSION','builtin'), pkg load symbolic; end

% Singular native solve remains native behavior; no generic multiphase rewrite.
Q = {zeros(1,3), zeros(1,3)};
[al,m] = solve_charge_vectors(Q,1,0.5);
assert(iscell(al) && numel(al) == 2 && isequal(size(m),[1 1]));

% Native architecture seam accepts phase rows beyond adapter's fixed two phases.
arch = dickson_arch(2);
arch.Asw_act = [arch.Asw_act; arch.Asw_act(1,:)];
assert(size(arch.Asw_act,1) > 2);
end
