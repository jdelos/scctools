#include "native_boundary.h"
#include <cassert>
#include <cstring>
#include <string>
int main() {
  char *ok = scctools_submit_json("{\"version\":1,\"architecture\":\"qfa-graph\",\"stages\":2,\"phases\":2,\"capacitors\":2,\"duty\":0.5}");
  assert(ok && std::strstr(ok,"\"type\":\"result\"") && std::strstr(ok,"\"A\"") && std::strstr(ok,"\"m\"")); scctools_free(ok);
  char *bad = scctools_submit_json("{\"version\":1,\"architecture\":\"qfa-graph\",\"stages\":2,\"phases\":3,\"capacitors\":2,\"duty\":0.5}");
  assert(bad && std::strstr(bad,"UNSUPPORTED_PHASE_COUNT")); scctools_free(bad);
  char *too_few = scctools_submit_json("{\"version\":1,\"architecture\":\"qfa-graph\",\"stages\":2,\"phases\":2,\"capacitors\":1,\"duty\":0.5}");
  assert(too_few && std::strstr(too_few,"INVALID_INPUT")); scctools_free(too_few);
  for (const char *invalid : {"{\"version\":1,\"architecture\":\"qfa-graph\",\"stages\":2,\"phases\":2,\"capacitors\":2}", " {\"version\":1,\"architecture\":\"qfa-graph\",\"stages\":2,\"phases\":2,\"capacitors\":2,\"duty\":0.5} trailing"}) {
    char *rejected = scctools_submit_json(invalid); assert(rejected && std::strstr(rejected,"INVALID_INPUT")); scctools_free(rejected);
  }
  return 0;
}
