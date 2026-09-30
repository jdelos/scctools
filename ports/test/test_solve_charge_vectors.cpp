#include "solve_charge_vectors.h"
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <symengine/parser.h>
#include <symengine/real_double.h>
#include <symengine/symbol.h>
using SymEngine::DenseMatrix;
static DenseMatrix matrix(unsigned rows, unsigned cols, const char *const *v) { DenseMatrix out(rows, cols); for (unsigned r=0;r<rows;++r) for(unsigned c=0;c<cols;++c) out.set(r,c,SymEngine::parse(v[r*cols+c])); return out; }
static DenseMatrix matrix(unsigned rows, unsigned cols, std::initializer_list<const char *> v) { assert(v.size() == rows * cols); DenseMatrix out(rows, cols); unsigned i = 0; for (const char *x : v) { out.set(i / cols, i % cols, SymEngine::parse(x)); ++i; } return out; }
static DenseMatrix duty(const SymEngine::RCP<const SymEngine::Symbol> &d) { DenseMatrix out(1,2); out.set(0,0,d); out.set(0,1,SymEngine::sub(SymEngine::integer(1),d)); return out; }
static void same(const DenseMatrix &a, const DenseMatrix &b) { assert(a.nrows()==b.nrows() && a.ncols()==b.ncols()); for (unsigned r=0;r<a.nrows();++r) for(unsigned c=0;c<a.ncols();++c) { if (!SymEngine::eq(*a.get(r,c),*b.get(r,c))) std::cerr << r << "," << c << "\n"; assert(SymEngine::eq(*a.get(r,c),*b.get(r,c))); } }
static void zero(const DenseMatrix &a) { for (unsigned r=0;r<a.nrows();++r) for (unsigned c=0;c<a.ncols();++c) assert(SymEngine::eq(*a.get(r,c), *SymEngine::integer(0))); }
int main() {
    auto d = SymEngine::symbol("D");
    const char *q1v[] = {"-1", "1", "2", "1", "0", "3"};
    const char *q2v[] = {"-1", "0", "1/2", "0", "1", "2"};
    bool rejected = false;
    try { (void)solve_charge_vectors({matrix(2,3,q1v), matrix(2,3,q2v)}, 1, duty(d), {d}); } catch (const std::invalid_argument &) { rejected = true; }
    assert(rejected);
    const char *q1short[] = {"-1", "1", "2"};
    auto result = solve_charge_vectors({matrix(1,3,q1short), matrix(2,3,q2v)}, 1, duty(d), {d});
    assert(result.Qx.nrows() == 4 && result.Qo.nrows() == 4);
    assert(result.a.size() == 2 && result.a[0].nrows() == 2 && result.a[1].nrows() == 2);
    assert(result.m.nrows() == 1);

    const char *multi0[] = {"-1", "1", "2", "3"};
    const char *multi1[] = {"-1", "0", "1/2", "2", "0", "1", "1", "4"};
    auto multi = solve_charge_vectors({matrix(1,4,multi0), matrix(2,4,multi1)}, 1, duty(d), {d});
    same(multi.Qx, matrix(4,4, {"1","1","0","0", "0","0","1","0", "0","0","0","1", "0","1","0","1"}));
    same(multi.a[0], matrix(2,2, {"-3","-7", "1","4"}));
    same(multi.a[1], matrix(2,2, {"-1/2","-2", "-1","-4"}));
    same(multi.m, matrix(2,1, {"-7/2","-9"}));
    for (unsigned p = 0; p < 2; ++p) {
        DenseMatrix residual(2, 2);
        for (unsigned r = 0; r < 2; ++r) for (unsigned c = 0; c < 2; ++c) {
            SymEngine::RCP<const SymEngine::Basic> value = SymEngine::integer(0);
            for (unsigned k = 0; k < 4; ++k) {
                const auto &phase = k < 2 ? multi.a[0] : multi.a[1];
                value = SymEngine::add(value, SymEngine::mul(multi.Qx.get(p * 2 + r, k), phase.get(k % 2, c)));
            }
            residual.set(r, c, SymEngine::add(value, multi.Qo.get(p * 2 + r, c)));
        }
        zero(residual);
    }
    const char *singular0[] = {"-1", "1", "0", "0"};
    const char *singular1[] = {"-1", "0", "1", "0", "-1", "0", "1", "0"};
    bool singular = false;
    try { (void)solve_charge_vectors({matrix(1,4,singular0), matrix(2,4,singular1)}, 1, duty(d), {d}); }
    catch (const std::runtime_error &) { singular = true; }
    assert(singular);

    DenseMatrix low(1,2); low.set(0,0,SymEngine::real_double(-0.1)); low.set(0,1,SymEngine::real_double(1.1));
    rejected = false; try { (void)solve_charge_vectors({matrix(1,3,q1v), matrix(2,3,q2v)}, 1, low, {}); } catch (const std::invalid_argument &) { rejected = true; } assert(rejected);
    DenseMatrix high(1,2); high.set(0,0,SymEngine::real_double(1.1)); high.set(0,1,SymEngine::real_double(-0.1));
    rejected = false; try { (void)solve_charge_vectors({matrix(1,3,q1v), matrix(2,3,q2v)}, 1, high, {}); } catch (const std::invalid_argument &) { rejected = true; } assert(rejected);
    rejected = false; try { (void)solve_charge_vectors({matrix(1,3,q1v), matrix(2,3,q2v)}, 1, duty(d), {SymEngine::RCP<const SymEngine::Symbol>()}); } catch (const std::invalid_argument &) { rejected = true; } assert(rejected);
}
