function test_characterize_dickson_two_phase()
reports = characterize_dickson_two_phase();
assert([reports.n_caps] == [2 3]);
for r = reports
    assert(all(r.Asw_act(1,1:2:end) == 1));
    assert(all(r.Asw_act(2,2:2:end) == 1));
    assert(isempty(r.A) == ~isempty(r.error));
    assert(isempty(r.m) == ~isempty(r.error));
end
end
