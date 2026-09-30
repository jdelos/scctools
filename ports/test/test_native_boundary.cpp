#include "native_boundary.h"
#include <cassert>
#include <cstring>
#include <string>
int main() {
  char *ok = scctools_submit_json(" { \"duty\": 0.5, \"capacitors\": 2, \"phases\": 2, \"stages\": 2, \"architecture\": \"qfa-graph\", \"version\": 1 } ");
  assert(ok && std::strstr(ok,"\"type\":\"result\"") && std::strstr(ok,"\"A\"") && std::strstr(ok,"\"m\"")); scctools_free(ok);
  char *bad = scctools_submit_json("{\"version\":1,\"architecture\":\"qfa-graph\",\"stages\":2,\"phases\":3,\"capacitors\":2,\"duty\":0.5}");
  assert(bad && std::strstr(bad,"UNSUPPORTED_PHASE_COUNT")); scctools_free(bad);
  char *too_few = scctools_submit_json("{\"version\":1,\"architecture\":\"qfa-graph\",\"stages\":2,\"phases\":2,\"capacitors\":1,\"duty\":0.5}");
  assert(too_few && std::strstr(too_few,"INVALID_INPUT")); scctools_free(too_few);
  for (const char *invalid : {"{\"version\":1,\"architecture\":\"qfa-graph\",\"stages\":2,\"phases\":2,\"capacitors\":2}", " {\"version\":1,\"architecture\":\"qfa-graph\",\"stages\":2,\"phases\":2,\"capacitors\":2,\"duty\":0.5} trailing"}) {
    char *rejected = scctools_submit_json(invalid); assert(rejected && std::strstr(rejected,"INVALID_INPUT")); scctools_free(rejected);
  }
  char *singular = scctools_submit_json("{\"version\":1,\"operation\":\"solve-charge-vectors\",\"capacitors\":1,\"duty\":[0.5,0.5],\"cutsets\":[[[0,0,0]],[[0,0,0],[0,0,0]]]}" );
  assert(singular && std::strstr(singular, "SINGULAR_SYSTEM")); scctools_free(singular);
  char *duplicate = scctools_submit_json("{\"version\":1,\"version\":1,\"architecture\":\"qfa-graph\",\"stages\":2,\"phases\":2,\"capacitors\":2,\"duty\":0.5}");
  assert(duplicate && std::strstr(duplicate, "INVALID_INPUT")); scctools_free(duplicate);
  return 0;
}
