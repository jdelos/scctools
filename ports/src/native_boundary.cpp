#include "native_boundary.h"
#include "dickson_hybrid_topology.h"
#include "solve_charge_vectors.h"
#include <symengine/symbol.h>
#include <symengine/real_double.h>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <stdexcept>
#include <cctype>

struct Json { enum Kind { OBJECT, ARRAY, STRING, NUMBER, BOOL, NIL } kind; std::map<std::string, Json> object; std::vector<Json> array; std::string string; double number{}; };
class JsonParser {
  const std::string &s; size_t p=0;
  void ws(){while(p<s.size() && std::isspace(static_cast<unsigned char>(s[p])))++p;}
  void expect(char c){ws(); if(p>=s.size()||s[p++]!=c) throw std::invalid_argument("malformed JSON");}
  std::string str(){ws(); expect('"'); std::string x; while(p<s.size()&&s[p]!='"'){ unsigned char c=static_cast<unsigned char>(s[p]); if(c<0x20) throw std::invalid_argument("control character in JSON string"); if(s[p]=='\\') throw std::invalid_argument("unsupported JSON escape"); x+=s[p++]; } expect('"'); return x;}
  Json value(){ws(); if(p>=s.size()) throw std::invalid_argument("malformed JSON"); if(s[p]=='{') return object(); if(s[p]=='[') return array(); if(s[p]=='"'){Json v{Json::STRING};v.string=str();return v;} if(s.compare(p,4,"true")==0){p+=4;return Json{Json::BOOL};} if(s.compare(p,5,"false")==0){p+=5;return Json{Json::BOOL};} if(s.compare(p,4,"null")==0){p+=4;return Json{Json::NIL};} char *e=nullptr; double n=std::strtod(s.c_str()+p,&e); if(e==s.c_str()+p) throw std::invalid_argument("malformed JSON"); p=e-s.c_str(); Json v{Json::NUMBER};v.number=n;return v; }
  Json object(){Json v{Json::OBJECT}; expect('{'); ws(); if(p<s.size()&&s[p]=='}'){++p;return v;} for(;;){std::string k=str(); if(v.object.count(k)) throw std::invalid_argument("duplicate member"); expect(':'); v.object.emplace(k,value()); ws(); if(p<s.size()&&s[p]=='}'){++p;return v;} expect(',');}}
  Json array(){Json v{Json::ARRAY}; expect('['); ws(); if(p<s.size()&&s[p]==']'){++p;return v;} for(;;){v.array.push_back(value());ws();if(p<s.size()&&s[p]==']'){++p;return v;}expect(',');}}
public: explicit JsonParser(const std::string &x):s(x){} Json parse(){Json v=value();ws();if(p!=s.size())throw std::invalid_argument("trailing JSON");return v;}
};
static char *reply(const std::string &s){char *p=static_cast<char*>(std::malloc(s.size()+1));if(p)std::memcpy(p,s.c_str(),s.size()+1);return p;}
static std::string error(const char *code,const std::string &message){return std::string("{\"version\":1,\"type\":\"error\",\"error\":{\"code\":\"")+code+"\",\"message\":\""+message+"\"}}";}
static const Json &field(const Json &o,const char *k){auto i=o.object.find(k);if(i==o.object.end())throw std::invalid_argument(std::string("missing field: ")+k);return i->second;}
static int integer(const Json &v,const char *name){if(v.kind!=Json::NUMBER||v.number!=static_cast<int>(v.number))throw std::invalid_argument(std::string("wrong type: ")+name);return static_cast<int>(v.number);}
static double number(const Json &v,const char *name){if(v.kind!=Json::NUMBER)throw std::invalid_argument(std::string("wrong type: ")+name);return v.number;}
static std::string matrix(const SymEngine::DenseMatrix &m){std::ostringstream out;out<<'[';for(unsigned r=0;r<m.nrows();++r){if(r)out<<',';out<<'[';for(unsigned c=0;c<m.ncols();++c){if(c)out<<',';out<<'"'<<*m.get(r,c)<<'"';}out<<']';}return out<<']',out.str();}
static SymEngine::DenseMatrix input_matrix(const Json &v){if(v.kind!=Json::ARRAY||v.array.empty())throw std::invalid_argument("cutset must be matrix");unsigned cols=0;for(const auto&r:v.array){if(r.kind!=Json::ARRAY||r.array.empty()||(cols&&r.array.size()!=cols))throw std::invalid_argument("cutset must be rectangular");cols=r.array.size();}SymEngine::DenseMatrix m(v.array.size(),cols);for(unsigned r=0;r<v.array.size();++r)for(unsigned c=0;c<cols;++c)m.set(r,c,SymEngine::real_double(number(v.array[r].array[c],"cutset entry")));return m;}
extern "C" char *scctools_submit_json(const char *request_json){try{if(!request_json)return reply(error("INVALID_INPUT","request is null"));Json req=JsonParser(request_json).parse();if(req.kind!=Json::OBJECT)throw std::invalid_argument("request must be object");for(const auto &x:req.object)if(x.first!="version"&&x.first!="architecture"&&x.first!="stages"&&x.first!="phases"&&x.first!="capacitors"&&x.first!="duty"&&x.first!="operation"&&x.first!="cutsets")throw std::invalid_argument("unknown field: "+x.first);if(integer(field(req,"version"),"version")!=1)throw std::invalid_argument("invalid version");
 if(req.object.count("operation")){const Json &op=field(req,"operation");if(op.kind!=Json::STRING||op.string!="solve-charge-vectors")throw std::invalid_argument("invalid operation");if(req.object.count("architecture")||req.object.count("stages")||req.object.count("phases"))throw std::invalid_argument("wrong operation fields");int caps=integer(field(req,"capacitors"),"capacitors");const Json &ds=field(req,"duty");if(ds.kind!=Json::ARRAY||ds.array.size()!=2)throw std::invalid_argument("duty must have two entries");SymEngine::DenseMatrix duty(1,2);for(unsigned i=0;i<2;++i)duty.set(0,i,SymEngine::real_double(number(ds.array[i],"duty")));const Json &cs=field(req,"cutsets");if(cs.kind!=Json::ARRAY||cs.array.size()!=2)throw std::invalid_argument("cutsets must have two phases");std::vector<SymEngine::DenseMatrix> q{input_matrix(cs.array[0]),input_matrix(cs.array[1])};std::vector<SymEngine::RCP<const SymEngine::Symbol>> symbols;ChargeSolution result=solve_charge_vectors(q,caps,duty,symbols);std::ostringstream out;out<<"{\"version\":1,\"type\":\"result\",\"operation\":\"solve-charge-vectors\",\"A\":["<<matrix(result.a[0])<<","<<matrix(result.a[1])<<"],\"m\":"<<matrix(result.m)<<"}";return reply(out.str());}
 const Json &arch=field(req,"architecture");if(arch.kind!=Json::STRING||arch.string!="qfa-graph")throw std::invalid_argument("invalid architecture");if(integer(field(req,"stages"),"stages")!=2)throw std::invalid_argument("invalid stages");if(integer(field(req,"phases"),"phases")!=2)return reply(error("UNSUPPORTED_PHASE_COUNT","only two phases supported"));int caps=integer(field(req,"capacitors"),"capacitors");if(caps<2||caps>3)throw std::invalid_argument("capacitors must be 2 or 3");double duty=number(field(req,"duty"),"duty");if(duty!=.25&&duty!=.5)throw std::invalid_argument("invalid duty");Topology t=dickson_hybrid_topology(caps,duty,true,false);std::ostringstream out;out<<"{\"version\":1,\"type\":\"result\",\"architecture\":\"qfa-graph\",\"parameters\":{\"stages\":2,\"phases\":2,\"capacitors\":"<<caps<<",\"duty\":"<<duty<<"},\"ordering\":{\"phases\":[0,1],\"capacitors\":"<<caps<<"},\"A\":[";for(unsigned p=0;p<t.phase.size();++p){if(p)out<<',';out<<matrix(t.phase[p].get_a_vector());}out<<"],\"m\":"<<matrix(t.m_ratios)<<",\"metadata\":{\"native\":true,\"provenance\":\"native-qfa\"}}";return reply(out.str());}catch(const std::exception&e){return reply(error(std::string(e.what()).find("singular")!=std::string::npos?"SINGULAR_SYSTEM":"INVALID_INPUT",e.what()));}}
extern "C" void scctools_free(char *p){std::free(p);}
