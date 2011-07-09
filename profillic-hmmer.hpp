#ifndef __GALOSH_PROFILLICHMMER_HPP__
#define __GALOSH_PROFILLICHMMER_HPP__

// HMMoC-BFloat-Algebra:
#include "Algebra.hpp"

// Galosh:
#include "Profile.hpp"
using galosh::ProfileTreeRoot;

/////////////// For profillic-hmmer //////////////////////////////////
/// Stuff we needed to modify in order to compile it in c++:
/// NOTE I had to change hmmer3/easel/esl_msa.h, where keyword "new" was being used as an argument name in a predeclaration for esl_msa_Copy (..).  It now reads:
//extern int      esl_msa_Copy (const ESL_MSA *msa, ESL_MSA *_new);
#define ESL_ALLOC_CPP(arg_type, p, size) do {            \
    if ( ( (p) = ( static_cast<arg_type *>( malloc(size)) ) ) == NULL) { \
       status = eslEMEM;\
       esl_exception(eslEMEM, __FILE__, __LINE__, "malloc of size %d failed", size);\
       goto ERROR;\
     }} while (0)
#define ESL_RALLOC_CPP(arg_type, p, tmp, newsize) do {           \
    if ((p) == NULL) { \
     (tmp) = ( malloc(newsize) );         \
    } else             { (tmp) = realloc((p), (newsize)); } \
     if ((tmp) != NULL) (p) = static_cast<arg_type *>(tmp);\
     else {\
       status = eslEMEM;\
       esl_exception(eslEMEM, __FILE__, __LINE__, "realloc for size %d failed", newsize);\
       goto ERROR;\
     }} while (0)
///
/////////////// End profillic-hmmer //////////////////////////////////

#endif // __GALOSH_PROFILLICHMMER_HPP__
