
#include <symengine/expression.h>
#include <symengine/symbol.h>
#include <string>
#include <emscripten/emscripten.h>

using SymEngine::Expression;
using SymEngine::symbol;

extern "C" {

// Function to compute the derivative of an expression
EMSCRIPTEN_KEEPALIVE
const char* differentiate(const char* expr_str, const char* var_str) {
    static std::string result_str;  // Persistent buffer

    try {
        // Define the variable and expression
        Expression var(symbol(var_str));
        Expression expr(expr_str);

        // Differentiate the expression
        Expression result = expr.diff(var);

        // Convert to string
        result_str = result.get_basic()->__str__();
    } catch (...) {
        result_str = "Error: Invalid expression.";
    }

    return result_str.c_str();
}

}
