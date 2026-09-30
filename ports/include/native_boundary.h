#ifndef NATIVE_BOUNDARY_H
#define NATIVE_BOUNDARY_H
#ifdef __cplusplus
extern "C" {
#endif
/* Returned strings belong to caller; release with scctools_free. */
char *scctools_submit_json(const char *request_json);
void scctools_free(char *response_json);
#ifdef __cplusplus
}
#endif
#endif
