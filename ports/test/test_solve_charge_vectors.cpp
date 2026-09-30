#include "solve_charge_vectors.h"
#include <cassert>
#include <stdexcept>
#include <symengine/parser.h>
#include <symengine/real_double.h>
#include <symengine/symbol.h>
using SymEngine::DenseMatrix;
static DenseMatrix matrix(unsigned rows, unsigned cols, const char *const *v) { DenseMatrix out(rows, cols); for (unsigned r=0;r<rows;++r) for(unsigned c=0;c<cols;++c) out.set(r,c,SymEngine::parse(v[r*cols+c])); return out; }
static DenseMatrix duty(const SymEngine::RCP<const SymEngine::Symbol> &d) { DenseMatrix out(1,2); out.set(0,0,d); out.set(0,1,SymEngine::sub(SymEngine::integer(1),d)); return out; }
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
    DenseMatrix low(1,2); low.set(0,0,SymEngine::real_double(-0.1)); low.set(0,1,SymEngine::real_double(1.1));
    rejected = false; try { (void)solve_charge_vectors({matrix(1,3,q1v), matrix(2,3,q2v)}, 1, low, {}); } catch (const std::invalid_argument &) { rejected = true; } assert(rejected);
    DenseMatrix high(1,2); high.set(0,0,SymEngine::real_double(1.1)); high.set(0,1,SymEngine::real_double(-0.1));
    rejected = false; try { (void)solve_charge_vectors({matrix(1,3,q1v), matrix(2,3,q2v)}, 1, high, {}); } catch (const std::invalid_argument &) { rejected = true; } assert(rejected);
    rejected = false; try { (void)solve_charge_vectors({matrix(1,3,q1v), matrix(2,3,q2v)}, 1, duty(d), {SymEngine::RCP<const SymEngine::Symbol>()}); } catch (const std::invalid_argument &) { rejected = true; } assert(rejected);
}
