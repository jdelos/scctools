#include "native_boundary.h"
#include "dickson_hybrid_topology.h"
#include <symengine/symbol.h>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <string>
#include <regex>

static char *reply(const std::string &s) { char *p = static_cast<char *>(std::malloc(s.size()+1)); if (p) std::memcpy(p,s.c_str(),s.size()+1); return p; }
static std::string error(const char *code, const std::string &message) {
    return std::string("{\"version\":1,\"type\":\"error\",\"error\":{\"code\":\"") + code + "\",\"message\":\"" + message + "\"}}";
}
static std::string matrix(const SymEngine::DenseMatrix &m) {
    std::ostringstream out; out << '[';
    for (unsigned r=0;r<m.nrows();++r) { if (r) out << ','; out << '['; for (unsigned c=0;c<m.ncols();++c) { if(c) out<<','; out << '\"' << *m.get(r,c) << '\"'; } out << ']'; }
    return out << ']', out.str();
}
extern "C" char *scctools_submit_json(const char *request_json) {
    try {
        if (!request_json) return reply(error("INVALID_INPUT","request is null"));
        const std::string req(request_json);
        std::smatch match;
        const std::regex schema(R"re(^\s*\{\s*"version"\s*:\s*(\d+)\s*,\s*"architecture"\s*:\s*"([^"]+)"\s*,\s*"stages"\s*:\s*(\d+)\s*,\s*"phases"\s*:\s*(\d+)\s*,\s*"capacitors"\s*:\s*(\d+)\s*,\s*"duty"\s*:\s*(0(?:\.25|\.5))\s*\}\s*$)re");
        if (!std::regex_match(req, match, schema)) return reply(error("INVALID_INPUT","request does not match version 1 schema"));
        if (match[1] != "1" || match[2] != "qfa-graph" || match[3] != "2") return reply(error("INVALID_INPUT","invalid architecture or stages"));
        if (match[4] != "2") return reply(error("UNSUPPORTED_PHASE_COUNT","only two phases supported"));
        int caps = std::stoi(match[5]);
        if (caps < 1 || caps > 3) return reply(error("INVALID_INPUT","capacitors must be 1, 2, or 3"));
        double duty = std::stod(match[6]);
        Topology t = dickson_hybrid_topology(caps, duty, true, false);
        std::ostringstream out;
        out << "{\"version\":1,\"type\":\"result\",\"architecture\":\"qfa-graph\",\"parameters\":{";
        out << "\"stages\":2,\"phases\":2,\"capacitors\":" << caps << ",\"duty\":" << duty << "},\"ordering\":{\"phases\":[0,1],\"capacitors\":" << caps << "},\"A\":[";
        for (unsigned p=0;p<t.phase.size();++p) { if(p) out<<','; out<<matrix(t.phase[p].get_a_vector()); }
        out << "],\"m\":" << matrix(t.m_ratios) << ",\"metadata\":{\"native\":true,\"provenance\":\"native-qfa\"}}";
        return reply(out.str());
    } catch (const std::exception &e) { return reply(error(std::string(e.what()).find("singular") != std::string::npos ? "SINGULAR_SYSTEM" : "INVALID_INPUT",e.what())); }
}
extern "C" void scctools_free(char *p) { std::free(p); }
