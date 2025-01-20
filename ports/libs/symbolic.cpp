#include <symengine/expression.h>
#include <symengine/symbol.h>
#include <iostream>

using namespace SymEngine;


extern "C" {
    const char* differentiate(const char* expr_str, const char* var_str) {
      static std::string result_str;     
      RCP<const Symbol> x = symbol(var_str);

      Expression expr(expr_str);
      Expression result = expr.diff(x);
      result_str = result.get_basic()->__str__();
      return result_str.c_str(); 
    }
}

