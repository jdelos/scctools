#include "native_boundary.h"
#include <iostream>
#include <string>
int main() {
  std::string request;
  while (std::getline(std::cin, request)) {
    char *result = scctools_submit_json(request.c_str());
    if (!result)
      return 1;
    std::cout << result << '\n';
    scctools_free(result);
  }
}
