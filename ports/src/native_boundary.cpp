#include "native_boundary.h"
#include "dickson_hybrid_topology.h"
#include <symengine/symbol.h>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <string>

static char *reply(const std::string &s) { char *p = static_cast<char *>(std::malloc(s.size()+1)); if (p) std::memcpy(p,s.c_str(),s.size()+1); return p; }
static std::string error(const char *code, const std::string &message) {
    return std::string("{\"version\":1,\"type\":\"error\",\"error\":{\"code\":\"") + code + "\",\"message\":\"" + message + "\"}}";
}
static std::string matrix(const SymEngine::DenseMatrix &m) {
    std::ostringstream out; out << '[';
    for (unsigned r=0;r<m.nrows();++r) { if (r) out << ','; out << '['; for (unsigned c=0;c<m.ncols();++c) { if(c) out<<','; out << '\"' << *m.get(r,c) << '\"'; } out << ']'; }
    return out << ']', out.str();
}
static bool has(const std::string &s, const char *needle) { return s.find(needle) != std::string::npos; }
extern "C" char *scctools_submit_json(const char *request_json) {
    try {
        if (!request_json) return reply(error("INVALID_INPUT","request is null"));
        std::string req(request_json);
        if (!has(req,"\"version\":1")) return reply(error("INVALID_INPUT","version must be 1"));
        if (!has(req,"\"architecture\":\"qfa-graph\"")) return reply(error("INVALID_INPUT","architecture must be qfa-graph"));
        if (!has(req,"\"phases\":2")) return reply(error("UNSUPPORTED_PHASE_COUNT","only two phases supported"));
        if (!has(req,"\"stages\":2")) return reply(error("INVALID_INPUT","stages must be 2"));
        int caps = has(req,"\"capacitors\":3") ? 3 : (has(req,"\"capacitors\":2") ? 2 : 0);
        if (!caps) return reply(error("INVALID_INPUT","capacitors must be 2 or 3"));
        bool numeric = has(req,"\"duty\":0.5") || has(req,"\"duty\":0.25");
        double duty = has(req,"\"duty\":0.25") ? .25 : .5;
        Topology t = dickson_hybrid_topology(caps, duty, true, false);
        std::ostringstream out;
        out << "{\"version\":1,\"type\":\"result\",\"architecture\":\"qfa-graph\",\"parameters\":{";
        out << "\"stages\":2,\"phases\":2,\"capacitors\":" << caps << ",\"duty\":" << duty << "},\"ordering\":{\"phases\":[0,1],\"capacitors\":" << caps << "},\"A\":[";
        for (unsigned p=0;p<t.phase.size();++p) { if(p) out<<','; out<<matrix(t.phase[p].get_a_vector()); }
        out << "],\"m\":" << matrix(t.m_ratios) << ",\"metadata\":{\"native\":true,\"numeric\":" << (numeric?"true":"false") << "}}";
        return reply(out.str());
    } catch (const std::exception &e) { return reply(error(has(e.what(),"singular") ? "SINGULAR_SYSTEM" : "INVALID_INPUT",e.what())); }
}
extern "C" void scctools_free(char *p) { std::free(p); }
