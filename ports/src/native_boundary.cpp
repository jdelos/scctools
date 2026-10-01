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
#include <cmath>
#include <limits>
#include <algorithm>
#include <symengine/integer.h>
#include <symengine/eval_double.h>
#include <symengine/visitor.h>

struct Json { enum Kind { OBJECT, ARRAY, STRING, NUMBER, BOOL, NIL } kind; std::map<std::string, Json> object; std::vector<Json> array; std::string string; double number{}; explicit Json(Kind k):kind(k){} };
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
static std::string error(const char *code,const std::string &message){std::string safe;for(unsigned char c:message){if(c=='"'||c=='\\')safe+='\\';if(c<0x20){safe+=' ';continue;}safe+=c;}return std::string("{\"version\":1,\"type\":\"error\",\"error\":{\"code\":\"")+code+"\",\"message\":\""+safe+"\"}}";}
static const Json &field(const Json &o,const char *k){auto i=o.object.find(k);if(i==o.object.end())throw std::invalid_argument(std::string("missing field: ")+k);return i->second;}
static int integer(const Json &v,const char *name){if(v.kind!=Json::NUMBER||!std::isfinite(v.number)||v.number<std::numeric_limits<int>::min()||v.number>std::numeric_limits<int>::max()||v.number!=std::trunc(v.number))throw std::invalid_argument(std::string("wrong type: ")+name);return static_cast<int>(v.number);}
static double number(const Json &v,const char *name){if(v.kind!=Json::NUMBER||!std::isfinite(v.number))throw std::invalid_argument(std::string("wrong type: ")+name);return v.number;}
static std::string matrix(const SymEngine::DenseMatrix &m){std::ostringstream out;out<<'[';for(unsigned r=0;r<m.nrows();++r){if(r)out<<',';out<<'[';for(unsigned c=0;c<m.ncols();++c){if(c)out<<',';out<<'"'<<*m.get(r,c)<<'"';}out<<']';}return out<<']',out.str();}
static SymEngine::DenseMatrix input_matrix(const Json &v){if(v.kind!=Json::ARRAY||v.array.empty())throw std::invalid_argument("cutset must be matrix");unsigned cols=0;for(const auto&r:v.array){if(r.kind!=Json::ARRAY||r.array.empty()||(cols&&r.array.size()!=cols))throw std::invalid_argument("cutset must be rectangular");cols=r.array.size();}SymEngine::DenseMatrix m(v.array.size(),cols);for(unsigned r=0;r<v.array.size();++r)for(unsigned c=0;c<cols;++c)m.set(r,c,SymEngine::real_double(number(v.array[r].array[c],"cutset entry")));return m;}
static std::string symbols(const std::vector<SymEngine::RCP<const SymEngine::Symbol>> &v) {
 std::ostringstream out; out << '['; for(unsigned i=0;i<v.size();++i){if(i)out<<',';out<<'"'<<*v[i]<<'"';} out<<']'; return out.str();
}
static std::string indexes(unsigned first,unsigned count) {
 std::ostringstream out; out<<'['; for(unsigned i=0;i<count;++i){if(i)out<<',';out<<first+i;} out<<']';return out.str();
}
static SymEngine::map_basic_basic substitutions(const Json &v,const Topology &t) {
 if(v.kind!=Json::OBJECT||v.object.size()!=5)throw std::invalid_argument("substitutions require five ordered vectors");
 SymEngine::map_basic_basic map;
 auto vector=[&](const char *name,const std::vector<SymEngine::RCP<const SymEngine::Symbol>> &syms,bool positive) {
  const auto &values=field(v,name); if(values.kind!=Json::ARRAY||values.array.size()!=syms.size())throw std::invalid_argument(std::string("wrong substitution length: ")+name);
  for(unsigned i=0;i<syms.size();++i){double x=number(values.array[i],name);if(positive?x<=0:x<0)throw std::invalid_argument(std::string("invalid substitution: ")+name);map[syms[i]]=SymEngine::real_double(x);}
 };
 const auto &d=field(v,"duty");if(d.kind!=Json::ARRAY||d.array.size()!=2)throw std::invalid_argument("duty substitution must have two entries");
 double d0=number(d.array[0],"duty"),d1=number(d.array[1],"duty");
 if(d0<=0||d1<=0||d0>=1||d1>=1||std::abs(d0+d1-1)>1e-12)throw std::invalid_argument("duty substitutions must be positive and sum to one");
 map[SymEngine::symbol("D")]=SymEngine::real_double(d0);
 vector("capacitances",t.capacitances,true);vector("switch_resistances",t.switch_resistances,false);vector("capacitor_esr",t.capacitor_esr,false);vector("frequency",{t.frequency},true);return map;
}
static std::string bundle(const Topology &t,const SymEngine::map_basic_basic &values={}) {
 auto emit=[&](const SymEngine::DenseMatrix &m){if(values.empty())return matrix(m);SymEngine::DenseMatrix n(m.nrows(),m.ncols());for(unsigned r=0;r<m.nrows();++r)for(unsigned c=0;c<m.ncols();++c){double x=SymEngine::eval_double(*m.get(r,c)->subs(values));if(!std::isfinite(x))throw std::invalid_argument("nonfinite substitution result");n.set(r,c,SymEngine::real_double(x));}return matrix(n);};
 std::ostringstream out;
 auto phases=[&](const char *name,SymEngine::DenseMatrix SCCPhase::*member){out<<'"'<<name<<"\":[";for(unsigned p=0;p<t.phase.size();++p){if(p)out<<',';out<<emit(t.phase[p].*member);}out<<"],";};
 phases("A",&SCCPhase::a_vector);phases("B",&SCCPhase::b_vector);phases("G",&SCCPhase::r_vector);phases("r",&SCCPhase::r_vector);phases("Ar",&SCCPhase::ar_vector);
 out<<"\"m\":"<<emit(t.m_ratios)<<",\"stress\":{\"capacitor_voltage\":"<<emit(t.vc)<<",\"switch_voltage_by_phase\":"<<emit(t.vr)<<",\"switch_current_by_output\":"<<emit(t.is)<<"},\"ZSSL\":"<<emit(t.ZSSL)<<",\"ZFSL\":"<<emit(t.ZFSL)<<",\"ZESR\":"<<emit(t.ZESR)<<",\"ZSCC\":"<<emit(t.ZSCC);return out.str();
}
extern "C" char *scctools_submit_json(const char *request_json){try{if(!request_json)return reply(error("INVALID_INPUT","request is null"));Json req=JsonParser(request_json).parse();if(req.kind!=Json::OBJECT)throw std::invalid_argument("request must be object");for(const auto &x:req.object)if(x.first!="version"&&x.first!="architecture"&&x.first!="stages"&&x.first!="phases"&&x.first!="capacitors"&&x.first!="duty"&&x.first!="operation"&&x.first!="cutsets"&&x.first!="substitutions")throw std::invalid_argument("unknown field: "+x.first);if(integer(field(req,"version"),"version")!=1)throw std::invalid_argument("invalid version");
 if(req.object.count("operation")){const Json &op=field(req,"operation");if(op.kind!=Json::STRING||op.string!="solve-charge-vectors")throw std::invalid_argument("invalid operation");if(req.object.count("architecture")||req.object.count("stages")||req.object.count("phases")||req.object.count("substitutions"))throw std::invalid_argument("wrong operation fields");int caps=integer(field(req,"capacitors"),"capacitors");const Json &ds=field(req,"duty");if(ds.kind!=Json::ARRAY||ds.array.size()!=2)throw std::invalid_argument("duty must have two entries");SymEngine::DenseMatrix duty(1,2);for(unsigned i=0;i<2;++i)duty.set(0,i,SymEngine::real_double(number(ds.array[i],"duty")));const Json &cs=field(req,"cutsets");if(cs.kind!=Json::ARRAY||cs.array.size()!=2)throw std::invalid_argument("cutsets must have two phases");std::vector<SymEngine::DenseMatrix> q{input_matrix(cs.array[0]),input_matrix(cs.array[1])};std::vector<SymEngine::RCP<const SymEngine::Symbol>> symbols;ChargeSolution result=solve_charge_vectors(q,caps,duty,symbols);std::ostringstream out;out<<"{\"version\":1,\"type\":\"result\",\"operation\":\"solve-charge-vectors\",\"A\":["<<matrix(result.a[0])<<","<<matrix(result.a[1])<<"],\"m\":"<<matrix(result.m)<<"}";return reply(out.str());}
 if(req.object.count("cutsets"))throw std::invalid_argument("wrong topology fields");
 const Json &arch=field(req,"architecture");if(arch.kind!=Json::STRING||arch.string!="qfa-graph")throw std::invalid_argument("invalid architecture");if(integer(field(req,"stages"),"stages")!=2)throw std::invalid_argument("invalid stages");if(integer(field(req,"phases"),"phases")!=2)return reply(error("UNSUPPORTED_PHASE_COUNT","only two phases supported"));int caps=integer(field(req,"capacitors"),"capacitors");if(caps<2||caps>3)throw std::invalid_argument("capacitors must be 2 or 3");double duty=number(field(req,"duty"),"duty");if(duty!=.25&&duty!=.5)throw std::invalid_argument("invalid duty");auto D=SymEngine::symbol("D");
 Topology t=dickson_hybrid_topology(caps,SymEngine::Expression(D),{D},true,false);
 std::ostringstream out;out<<"{\"version\":1,\"type\":\"result\",\"architecture\":\"qfa-graph\",\"parameters\":{\"stages\":2,\"phases\":2,\"capacitors\":"<<caps<<",\"duty\":"<<duty<<"},\"ordering\":{\"index_base\":1,\"phases\":[1,2],\"nodes\":"<<indexes(1,t.m_ratios.nrows()+1)<<",\"outputs\":"<<indexes(2,t.m_ratios.nrows())<<",\"capacitors\":"<<indexes(1,caps)<<",\"switches\":"<<indexes(1,t.N_sw)<<",\"A_rows\":[\"supply\",\"capacitors\"],\"B_rows\":\"capacitors\",\"r_rows\":\"capacitors\",\"Ar_rows\":\"capacitors then active switches\",\"active_switches\":[";
 for(unsigned p=0;p<2;++p){if(p)out<<',';out<<'[';for(unsigned i=0;i<t.phase[p].sw_idxs.size();++i){if(i)out<<',';out<<t.phase[p].sw_idxs[i];}out<<']';}
 auto ordered=t.ordered_symbols;ordered.insert(ordered.end(),t.capacitances.begin(),t.capacitances.end());ordered.insert(ordered.end(),t.switch_resistances.begin(),t.switch_resistances.end());ordered.insert(ordered.end(),t.capacitor_esr.begin(),t.capacitor_esr.end());ordered.push_back(t.frequency);
 std::sort(ordered.begin(),ordered.end(),[](const auto &a,const auto &b){return a->get_name()<b->get_name();});
 out<<"]},\"symbols\":{\"duty\":[\"D\",\"1-D\"],\"capacitances\":"<<symbols(t.capacitances)<<",\"switch_resistances\":"<<symbols(t.switch_resistances)<<",\"capacitor_esr\":"<<symbols(t.capacitor_esr)<<",\"frequency\":[\"fsw\"],\"matlab_compatible\":"<<symbols(ordered)<<"},"<<bundle(t);
 if(req.object.count("substitutions")) {try{out<<",\"evaluated\":{"<<bundle(t,substitutions(field(req,"substitutions"),t))<<'}';}catch(const std::exception &e){return reply(error("INVALID_SUBSTITUTION",e.what()));}}
 out<<",\"metadata\":{\"native\":true,\"provenance\":\"native-qfa\",\"ZSCC\":\"elementwise analytical root-sum-square approximation; not transient impedance\",\"ZFSL\":\"includes capacitor ESR; ZESR is its capacitor-only contribution\"}}";return reply(out.str());}catch(const std::exception&e){return reply(error(std::string(e.what()).find("singular")!=std::string::npos?"SINGULAR_SYSTEM":"INVALID_INPUT",e.what()));}}
extern "C" void scctools_free(char *p){std::free(p);}
