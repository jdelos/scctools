#include "solve_charge_vectors.h"
#include <cassert>
#include <stdexcept>
#include <symengine/parser.h>
#include <symengine/symbol.h>

using SymEngine::DenseMatrix;
using SymEngine::Expression;

static DenseMatrix matrix(unsigned rows, unsigned cols, const char *const *values) {
    DenseMatrix out(rows, cols);
    for (unsigned r = 0; r < rows; ++r)
        for (unsigned c = 0; c < cols; ++c)
            out.set(r, c, SymEngine::parse(values[r * cols + c]));
    return out;
}
static void expect(const SymEngine::Basic &actual, const char *value) {
    auto expected = SymEngine::parse(value);
    assert(SymEngine::eq(actual, *expected));
}
static void expect_matrix(const DenseMatrix &actual, unsigned rows, unsigned cols,
                          const char *const *values) {
    assert(actual.nrows() == rows && actual.ncols() == cols);
    for (unsigned r = 0; r < rows; ++r)
        for (unsigned c = 0; c < cols; ++c)
            expect(*actual.get(r, c), values[r * cols + c]);
}
static DenseMatrix duty(const SymEngine::RCP<const SymEngine::Symbol> &d) {
    DenseMatrix out(1, 2);
    out.set(0, 0, d);
    out.set(0, 1, SymEngine::sub(SymEngine::integer(1), d));
    return out;
}
static DenseMatrix fixture_q1() {
    const char *v[] = {"-1", "1", "2", "1", "0", "3"};
    return matrix(2, 3, v);
}
static DenseMatrix fixture_q2() {
    const char *v[] = {"-1", "0", "1/2", "0", "1", "2"};
    return matrix(2, 3, v);
}
int main() {
    auto d = SymEngine::symbol("D");
    auto result = solve_charge_vectors({fixture_q1(), fixture_q2()}, 1, duty(d), {d});
    const char *qx[] = {"1", "1", "0", "0", "-1", "0", "0", "0", "0", "0", "1", "0", "0", "1", "0", "1"};
    const char *qo[] = {"2", "3", "1/2", "2"};
    expect_matrix(result.Qx, 4, 4, qx);
    expect_matrix(result.Qo, 4, 1, qo);
    assert(result.a.size() == 2);
    const char *a1[] = {"3", "-5"};
    const char *a2[] = {"-1/2", "3"};
    expect_matrix(result.a[0], 2, 1, a1);
    expect_matrix(result.a[1], 2, 1, a2);
    expect(*result.m.get(0, 0), "5/2");

    auto symbolic = solve_charge_vectors({fixture_q1(), matrix(2, 3, (const char *[]){"-1", "0", "D", "0", "1", "4*D"})}, 1, duty(d), {d});
    const char *sa1[] = {"3", "-5"};
    const char *sa2[] = {"-D", "5-4*D"};
    expect_matrix(symbolic.a[0], 2, 1, sa1);
    expect_matrix(symbolic.a[1], 2, 1, sa2);
    expect(*symbolic.m.get(0, 0), "3-D");

    bool rejected = false;
    try { (void)solve_charge_vectors({fixture_q1()}, 1, duty(d), {d}); } catch (const std::invalid_argument &) { rejected = true; }
    assert(rejected);
    DenseMatrix bad_duty(1, 2); bad_duty.set(0, 0, d); bad_duty.set(0, 1, d);
    rejected = false;
    try { (void)solve_charge_vectors({fixture_q1(), fixture_q2()}, 1, bad_duty, {d}); } catch (const std::invalid_argument &) { rejected = true; }
    assert(rejected);
    rejected = false;
    try { (void)solve_charge_vectors({fixture_q1(), fixture_q2()}, 1, duty(d), {}); } catch (const std::invalid_argument &) { rejected = true; }
    assert(rejected);

    DenseMatrix null_duty(1, 2); null_duty.set(0, 0, d);
    rejected = false;
    try { (void)solve_charge_vectors({fixture_q1(), fixture_q2()}, 1, null_duty, {d}); } catch (const std::invalid_argument &) { rejected = true; }
    assert(rejected);

    const char *unequal_q1_values[] = {"-1", "1", "1", "0", "1", "1"};
    const char *unequal_q2_values[] = {"-1", "0", "1"};
    auto unequal = solve_charge_vectors({matrix(2, 3, unequal_q1_values), matrix(1, 3, unequal_q2_values)}, 1, duty(d), {d});
    expect_matrix(unequal.Qx, 4, 4, (const char *[]) {"1", "1", "0", "0", "0", "1", "0", "0", "0", "0", "1", "0", "0", "1", "0", "1"});
    expect_matrix(unequal.Qo, 4, 1, (const char *[]) {"1", "1", "1", "0"});

    rejected = false;
    auto zero_q = matrix(2, 3, (const char *[]) {"0", "0", "0", "0", "0", "0"});
    try { (void)solve_charge_vectors({fixture_q1(), zero_q}, 1, duty(d), {d}); } catch (const std::runtime_error &) { rejected = true; }
    assert(rejected);
}
