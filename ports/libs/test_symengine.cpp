#include <symengine/expression.h>
#include <symengine/symbol.h>
#include <iostream>

using SymEngine::Expression;
using SymEngine::symbol;

int main() {
    Expression x(symbol("x"));
    Expression expr = x * x + 3 * x + 2;

    std::cout << "Expression: " << expr << std::endl;
    std::cout << "Derivative: " << expr.diff(x) << std::endl;

    return 0;
}

