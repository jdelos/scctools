#include "native_boundary.h"
#include <cassert>
#include <cstring>
#include <string>
int main() {
  char *ok = scctools_submit_json(" { \"duty\": 0.5, \"capacitors\": 2, \"phases\": 2, \"stages\": 2, \"architecture\": \"qfa-graph\", \"version\": 1 } ");
  assert(ok && std::strstr(ok,"\"type\":\"result\"") && std::strstr(ok,"\"A\"") && std::strstr(ok,"\"m\"")); for (const char *key : {"\"B\"", "\"r\"", "\"Ar\"", "\"ZSSL\"", "\"ZFSL\"", "\"ZESR\"", "\"ZSCC\"", "\"stress\"", "\"symbols\""}) assert(std::strstr(ok, key)); scctools_free(ok);
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
  char *control = scctools_submit_json("{\"bad\nfield\":1}");
  assert(control && std::strstr(control, "INVALID_INPUT") && !std::strchr(control, '\n')); scctools_free(control);
  const std::string prefix="{\"version\":1,\"architecture\":\"qfa-graph\",\"stages\":2,\"phases\":2,\"capacitors\":2,\"duty\":0.5,\"substitutions\":";
  for (const char *values : {"{}", "null", "[]", "{\"duty\":[0,1],\"capacitances\":[1,2],\"switch_resistances\":[1,1,1,1],\"capacitor_esr\":[0,0],\"frequency\":[1]}"}) {
    char *result=scctools_submit_json((prefix+values+"}").c_str());
    assert(result && std::strstr(result,"INVALID_SUBSTITUTION"));scctools_free(result);
  }
  char *huge=scctools_submit_json("{\"version\":1e300}");assert(huge && std::strstr(huge,"INVALID_INPUT"));scctools_free(huge);
  return 0;
}
