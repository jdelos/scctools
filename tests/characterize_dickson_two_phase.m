function reports = characterize_dickson_two_phase()
%CHARACTERIZE_DICKSON_TWO_PHASE Capture ordered two-phase Dickson inputs/outputs.
% No numeric A or m is substituted when QFA generation fails.
if exist('OCTAVE_VERSION','builtin')
    pkg load symbolic
end
reports = repmat(struct('n_caps', [], 'Acaps', [], 'Asw', [], ...
    'Asw_act', [], 'duties', [], 'A', [], 'm', [], 'error', '', ...
    'provenance', '', 'units', '', 'assumptions', '', 'tolerance', ''), 1, 2);
for k = 1:2
    n_caps = k + 1;
    arch = dickson_arch(n_caps);
    r = reports(k);
    r.n_caps = n_caps;
    r.Acaps = arch.Acaps;
    r.Asw = arch.Asw;
    r.Asw_act = arch.Asw_act;
    r.duties = {'D', '1-D'};
    r.A = [];
    r.m = [];
    r.provenance = 'dickson_arch.m and generic_switched_capacitor_class.m';
    r.units = 'incidence matrices dimensionless; duties fractions; A and m normalized';
    r.assumptions = 'ordered columns interleave phase 1, phase 2; symbolic D; exact algebra';
    r.tolerance = 'exact integer matrix comparison; no numeric tolerance';
    try
        model = generic_switched_capacitor_class(arch, 'Mode', 1, 'Duty', sym('D'));
        r.A = model.a;
        r.m = model.m;
    catch err
        r.error = sprintf('%s: %s', err.identifier, err.message);
    end
    reports(k) = r;
    fprintf('n_caps=%d\n', n_caps);
    fprintf('Acaps=\n'); disp(r.Acaps);
    fprintf('Asw=\n'); disp(r.Asw);
    fprintf('Asw_act=\n'); disp(r.Asw_act);
    fprintf('duties=[D, 1-D]\n');
    if isempty(r.error)
        fprintf('A=\n'); disp(r.A);
        fprintf('m=\n'); disp(r.m);
    else
        fprintf('A=<not generated>; m=<not generated>\nerror=%s\n', r.error);
    end
end
end
