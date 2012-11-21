/**
 * \file profillic-esl_msafile.hpp
 * \brief
 * Multiple sequence alignment file i/o (for profillic profiles)
 * \details
 * <pre>
 * Table of contents:
 *     1. Opening/closing an ESLX_MSAFILE.
 *     2. ESLX_MSAFILE_FMTDATA: optional added constraints on formats.
 *     3. Guessing file formats.
 *     4. Guessing alphabets.
 *     5. Random MSA flatfile access. [augmentation: ssi]
 *     6. Reading an MSA from an ESLX_MSAFILE.
 *     7. Writing an MSA to a stream.
 *     8. Utilities used by specific format parsers.
 *     9. Unit tests.
 *    10. Test driver.
 *    11. Examples.
 *    12. Copyright and license.
 * </pre>
 */
#ifndef __GALOSH_PROFILLICESLMSAFILE_HPP__
#define __GALOSH_PROFILLICESLMSAFILE_HPP__

extern "C" {
#include "esl_config.h"
}

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

extern "C" {
#include "easel.h"
#include "esl_mem.h"
#include "esl_msafile.h"
#ifdef eslAUGMENT_ALPHABET
#include "esl_alphabet.h"
#endif
}

/* /////////////// For profillic-hmmer ////////////////////////////////// */
/// \note TAH 8/12 Avoid the C++ keyword "new" in esl_msa.h
#define new _new
#include "profillic-hmmer.hpp"
extern "C" {
#include "esl_msa.h"
}
#undef new
#define eslMSAFILE_PROFILLIC       98103  /* A galosh profile (from profillic)   */
#define PRId64 "d"

// Predecs
int
profillic_eslx_msafile_Read(ESLX_MSAFILE *afp, ESL_MSA **ret_msa);

template <typename ProfileType>
int
profillic_eslx_msafile_Read(ESLX_MSAFILE *afp, ESL_MSA **ret_msa, ProfileType * profile_ptr );

template <typename ProfileType>
static int
profillic_esl_msafile_profile_Read(ESLX_MSAFILE *afp, ESL_MSA **ret_msa, ProfileType * profile_ptr );

/* /////////////// End profillic-hmmer ////////////////////////////////// */


/*****************************************************************
 *# 1. Opening/closing an ESLX_MSAFILE
 *****************************************************************/

static int profillic_msafile_Create    (ESLX_MSAFILE **ret_afp);
static int profillic_msafile_OpenBuffer(ESL_ALPHABET **byp_abc, ESL_BUFFER *bf, int format, ESLX_MSAFILE_FMTDATA *fmtd, ESLX_MSAFILE *afp);

/**
 * <pre>
 * Function:  eslx_msafile_Open()
 * Synopsis:  Open a multiple sequence alignment file for input.
 *
 * Purpose:   Open a multiple sequence alignment file <msafile> for input.
 *            Return an open <ESLX_MSAFILE> handle in <*ret_afp>.
 *
 *            <msafile> is usually the name of a file. Alignments may
 *            also be read from standard input, or from
 *            gzip-compressed files.  If <msafile> is ``-'', alignment
 *            input is taken from the standard input stream. If
 *            <msafile> ends in ``.gz'', alignment input is read
 *            through a pipe from <gzip -dc>.
 *            
 *            <byp_abc>, <env>, <format>, and <fmtd> support a variety
 *            of optional/advanced operations, as described
 *            below. Minimally, a caller can set <byp_abc> to <NULL>,
 *            <format> to <eslMSAFILE_UNKNOWN>, and <fmtd> to <NULL>,
 *            and <msafile> will be opened in text mode; in the
 *            current working directory; and its format will be
 *            autodetected.
 *            
 *            The <byp_abc> argument controls whether data are to be
 *            read in text or digital mode. In digital mode, alignment
 *            data are immediately digitized into an Easel internal
 *            alphabet (which among other things, allows various
 *            things to operate on sequence data more efficiently) and
 *            because an expected alphabet is known, parsers are able
 *            to detect invalid characters. The caller may either
 *            provide an alphabet (thus asserting what it's expected
 *            to be), or have <eslx_msafile_Open()> look at the file
 *            and guess what alphabet it appears to be (DNA or amino
 *            acid code, usually).  In text mode, alignment data are
 *            read verbatim. It might be advantageous for an
 *            application to read in text mode -- for example, if a
 *            variant alignment format is using characters in some
 *            special way, and you need to deal with them specially.
 *            All this goes through the setting of the passed-by-reference
 *            alphabet pointer <byp_abc>. If caller passes NULL for
 *            the <byp_abc> argument, input is in text mode. If caller
 *            provides a valid non-NULL <byp_abc> pointer but
 *            <*byp_abc> is NULL (that is, caller has declared
 *            <ESL_ALPHABET *abc = NULL> and passed <&abc> as an
 *            argument), then we attempt to guess the digital alphabet
 *            using <eslx_msafile_GuessAlphabet()>, based on the first
 *            alignment in the input. In this case, the new alphabet
 *            is allocated here and returned to the caller. If caller
 *            provides a digital alphabet (that is, <ESL_ALPHABET *abc
 *            = esl_alphabet_Create...()> and passed <&abc>), that's
 *            the alphabet we use.
 * 
 *            The <env> argument controls where we search for the
 *            <msafile>.  If <env> is <NULL>, only the current working
 *            directory is checked.  Optionally, caller can provide in
 *            <env> the name of an environment variable ("PFAMDB",
 *            perhaps), in which the routine can find a
 *            colon-delimited list of directories.  Then, if <msafile>
 *            is not found in the current working directory, we look
 *            for it in these directories, in the order they're
 *            listed.
 *
 *            The <format> argument allows the caller to either allow
 *            <eslx_msafile_Open()> to autodetect the file format of
 *            <msafile>, or to assert that it knows the file is in a
 *            particular format. If <format> is <eslMSAFILE_UNKNOWN>,
 *            format autodetection is performed. Other valid codes include:
 *             | <eslMSAFILE_STOCKHOLM>   | Stockholm format                    |
 *             | <eslMSAFILE_AFA>         | Aligned FASTA format                | 
 *             | <eslMSAFILE_CLUSTAL>     | Clustal format (strict)             |
 *             | <eslMSAFILE_CLUSTALLIKE> | Clustal-like  (MUSCLE, PROBCONS...) |
 *             | <eslMSAFILE_PHYLIP>      | PHYLIP interleaved format           |
 *             | <eslMSAFILE_PHYLIPS>     | PHYLIP sequential format            |
 *             | <eslMSAFILE_A2M>         | UCSC SAM A2M (dotless or dotful)    |
 *             | <eslMSAFILE_PSIBLAST>    | NCBI PSI-BLAST                      |
 *             | <eslMSAFILE_SELEX>       | a general alignment block format    |
 *
 *            The <fmtd> argument is an optional pointer to a
 *            <ESLX_MSAFILE_FMTDATA> structure that the caller may
 *            initialize and provide, in order to assert any
 *            additional unusual constraints on the input format --
 *            for example, to dictate that a PHYLIP format file has
 *            some nonstandard name field width. Generally, though,
 *            <fmtd> will be <NULL>.
 *
 * Args:      byp_abc   - digital alphabet to use, or NULL for text mode
 *                        if <*byp_abc> is NULL, guess the digital alphabet,
 *                        create it, and return it in <*byp_abc>.
 *                        If <*byp_abc> is a digital alphabet, use it.
 *            msafile   - name of alignment input to open;
 *                        if "-", read standard input;
 *                        if "*.gz", read through a <gzip -dc> pipe.
 *            env       - <NULL>, or the name of an environment variable
 *                        containing colon-delimited list of directories
 *                        in which to search for <msafile> (e.g. "PFAMDB").
 *            format    - format code, such as <eslMSAFILE_STOCKHOLM>;
 *                        or <eslMSAFILE_UNKNOWN> to autodetect format.
 *            fmtd      - <NULL>, or a pointer to an initialized 
 *                        <ESLX_MSAFILE_FMTDATA> structure, containing
 *                        any additional unusual constraints to apply
 *                        to the input format.
 *            *ret_afp  - RETURN: open MSA input stream.
 *
 * Returns:   <eslOK> on success, and <*ret_afp> is the newly opened msa file.
 *
 *            <eslENOTFOUND> if <msafile> doesn't exist or can't be
 *            opened for reading; or (in the case of a <.gz> file) if
 *            a <gzip> executable doesn't exist in user's <PATH> or
 *            can't be executed. <afp->errmsg> is something like 
 *            "couldn't open %s for reading", with <%s> being the 
 *            name of the msafile.
 *
 *            <eslENOFORMAT> if we tried to autodetect the file format
 *            (caller provided <format=eslMSAFILE_UNKNOWN>), and
 *            failed. <afp->errmsg> is something like "couldn't
 *            determine alignment input format".
 *
 *            <eslENOALPHABET> if we tried to autodetect the alphabet
 *            (caller provided <&abc>, <abc=NULL> to request digital
 *            mode w/ alphabet autodetection) but the alphabet could
 *            not be reliably guessed.
 *            
 *            <eslFAIL> in the case of a <.gz> file and the <gzip -dc>
 *            command fails on it.
 *            
 *            On any of these normal errors, <*ret_afp> is returned in
 *            an error state, containing a user-directed error message
 *            in <afp->errmsg> and (if relevant) the full path to
 *            <msafile> that we attempted to open in
 *            <afp->bf->filename>. See <eslx_msafile_OpenFailure()> for
 *            a function that gives a standard way of reporting these
 *            diagnostics to <stderr>.
 *            
 * Throws:    <eslEMEM> on allocation failure.
 *            <eslESYS> on a system call failure, such as <fread()>.
 *            <eslEINVAL> if we tried to use <stdin> but the <stdin> stream was
 *            invalid (in an error state, <NULL>, or at <EOF>).
 *            On thrown exceptions, <*ret_afp> is <NULL>.
 * </pre>
 */
int
profillic_eslx_msafile_Open(ESL_ALPHABET **byp_abc, const char *msafile, const char *env, int format, ESLX_MSAFILE_FMTDATA *fmtd, ESLX_MSAFILE **ret_afp)
{
  ESLX_MSAFILE *afp = NULL;
  int           status;

  if ( (status = profillic_msafile_Create(&afp)) != eslOK) goto ERROR;

  if ((status = esl_buffer_Open(msafile, env, &(afp->bf))) != eslOK)
    ESL_XFAIL(status, afp->errmsg, "%s", afp->bf->errmsg); /* ENOTFOUND; FAIL are normal here */

  if ( (status = profillic_msafile_OpenBuffer(byp_abc, afp->bf, format, fmtd, afp)) != eslOK) goto ERROR;

  *ret_afp = afp; 
  return eslOK;

 ERROR:  /* on normal errors, afp is returned in an error state */
  if (status == eslENOTFOUND || status == eslFAIL || status == eslEFORMAT || status == eslENODATA || eslENOALPHABET) 
    { afp->abc = NULL; *ret_afp = afp;}
  else 
    { if (afp) eslx_msafile_Close(afp);  *ret_afp = NULL; }
  return status;
}

static int
profillic_msafile_Create(ESLX_MSAFILE **ret_afp)
{
  ESLX_MSAFILE *afp = NULL;
  int           status;

  ESL_ALLOC_CPP(ESLX_MSAFILE, afp, sizeof(ESLX_MSAFILE));
  afp->bf         = NULL;
  afp->line       = NULL;
  afp->n          = 0;
  afp->linenumber = 0;
  afp->lineoffset = 0;
  afp->format     = eslMSAFILE_UNKNOWN;
  afp->abc        = NULL;
  afp->ssi        = NULL;
  afp->errmsg[0]  = '\0';

  eslx_msafile_fmtdata_Init(&(afp->fmtd));

  *ret_afp = afp;
  return eslOK;

 ERROR:
  *ret_afp = NULL;
  return status;
}

/* All input sources funnel through here.
 * Here, <afp> is already allocated and initialized, and the input
 * <bf> is opened successfully.
 */
static int
profillic_msafile_OpenBuffer(ESL_ALPHABET **byp_abc, ESL_BUFFER *bf, int format, ESLX_MSAFILE_FMTDATA *fmtd,  ESLX_MSAFILE *afp)
{
  ESL_ALPHABET        *abc       = NULL;
  int                  alphatype = eslUNKNOWN;
  int                  status;

  /* if caller provided <fmtd>, copy it into afp->fmtd */
  if (fmtd) eslx_msafile_fmtdata_Copy(fmtd, &(afp->fmtd));

  /* Determine the format */
  if (format == eslMSAFILE_UNKNOWN) 
    {
      status = eslx_msafile_GuessFileFormat(afp->bf, &format, &(afp->fmtd));
      if      (status == eslENOFORMAT) ESL_XFAIL(eslENOFORMAT, afp->errmsg, "couldn't determine alignment input format"); /* ENOFORMAT is normal failure */
      else if (status != eslOK)        goto ERROR;
    }
  afp->format = format;

  /* Determine the alphabet; set <abc>. (<abc> == NULL means text mode.)  */
  /* Note that GuessAlphabet() functions aren't allowed to use the inmap, because it isn't set yet */
#ifdef eslAUGMENT_ALPHABET
  if (byp_abc && *byp_abc)	/* Digital mode, and caller provided the alphabet */
    { 
      abc       = *byp_abc;
      alphatype = abc->type;
    } 
  else if (byp_abc)		/* Digital mode, and caller wants us to guess and create an alphabet */
    {
      status = eslx_msafile_GuessAlphabet(afp, &alphatype);
      if      (status == eslENOALPHABET) ESL_XFAIL(eslENOALPHABET, afp->errmsg, "couldn't guess alphabet (maybe try --dna/--rna/--amino if available)");
      else if (status != eslOK)          goto ERROR;
      if ( (abc = esl_alphabet_Create(alphatype))                == NULL) { status = eslEMEM; goto ERROR; }
    }    
#endif
  if (abc && ! byp_abc) ESL_EXCEPTION(eslEINCONCEIVABLE, "Your version of Easel does not include digital alphabet code."); 
  /* ^^^^^^^^^^^^^^^^^  this test interacts tricksily with the #ifdef above */
  afp->abc = abc;	/* with afp->abc set, the inmap config functions know whether to do digital/text    */

  /**
   * <pre>
   * Configure the format-specific, digital or text mode character
   * input map in afp->inmap.
   * All of these must:
   *    
   *    set inmap[0] to an appropriate 'unknown' character, to replace
   *       invalid input with.
   *    set ' ' to eslDSQ_IGNORE (if we're supposed to accept and skip
   *       it), or map it to a gap, or set it as eslDSQ_ILLEGAL.
   *    in digital mode, copy the abc->inmap
   *    in text mode, decide if we should accept most any
   *        non-whitespace character (isgraph()), or if the format is
   *        inherently restrictive and we should go with isalpha() +
   *        some other valid characters "_-.~*" instead.
   * </pre>
   */
  switch (afp->format) {
  case eslMSAFILE_A2M:          status = esl_msafile_a2m_SetInmap(      afp); break;
  case eslMSAFILE_AFA:          status = esl_msafile_afa_SetInmap(      afp); break;
  case eslMSAFILE_CLUSTAL:      status = esl_msafile_clustal_SetInmap(  afp); break;
  case eslMSAFILE_CLUSTALLIKE:  status = esl_msafile_clustal_SetInmap(  afp); break;
  case eslMSAFILE_PFAM:         status = esl_msafile_stockholm_SetInmap(afp); break;
  case eslMSAFILE_PHYLIP:       status = esl_msafile_phylip_SetInmap(   afp); break;
  case eslMSAFILE_PHYLIPS:      status = esl_msafile_phylip_SetInmap(   afp); break;
  case eslMSAFILE_PSIBLAST:     status = esl_msafile_psiblast_SetInmap( afp); break;
  case eslMSAFILE_SELEX:        status = esl_msafile_selex_SetInmap(    afp); break;
  case eslMSAFILE_STOCKHOLM:    status = esl_msafile_stockholm_SetInmap(afp); break;
  case eslMSAFILE_PROFILLIC:    status = eslOK;                               break; /// \todo status = profillic_esl_msafile_profile_SetInmap(afp); */ break;
  default: ESL_XEXCEPTION(eslENOFORMAT, "no such alignment file format");     break;
  }

  if (esl_byp_IsReturned(byp_abc)) *byp_abc = abc;
  return eslOK;

 ERROR:  /* on normal errors, afp is returned in an error state */
  if (abc && ! esl_byp_IsProvided(byp_abc)) { esl_alphabet_Destroy(abc); }
  if (esl_byp_IsReturned(byp_abc)) *byp_abc = NULL;
  afp->abc = NULL;
  return status;
}
/*------------- end, open/close an ESLX_MSAFILE -----------------*/

/*****************************************************************
 *# 6. Reading MSAs from input
 *****************************************************************/

/**
 * <pre>
 * Function:  eslx_msafile_Read()
 * Synopsis:  Read next MSA from input.
 *
 * Purpose:   Reads the next MSA from open MSA input <afp>, and return it in 
 *            <*ret_msa>.
 *
 * Args:      afp      - open alignment input stream
 *            *ret_msa - RETURN: alignment
 *
 * Returns:   <eslOK> on success. 
 *
 *            <eslEFORMAT> on a parse error, and <afp->errmsg> is set
 *            to a user-directed error message; <*ret_msa> is <NULL>.
 *
 *            If no alignment is found at all, returns <eslEOF>,
 *            and <afp->errmsg> is blank; <*ret_msa> is <NULL>.
 *
 *            On normal error, <afp> and the return status code may be
 *            passed to <eslx_msafile_ReadFailure()> to print diagnostics
 *            to <stderr> (including input source information and line
 *            number) and exit.
 *
 * Throws:    <eslEMEM> - an allocation failed.
 *            <eslESYS> - a system call such as fread() failed
 *            <eslEINCONCEIVABLE> - "impossible" corruption
 * </pre> 
 */
int
profillic_eslx_msafile_Read(ESLX_MSAFILE *afp, ESL_MSA **ret_msa)
{
  return profillic_eslx_msafile_Read(afp, ret_msa, (galosh::ProfileTreeRoot<seqan::Dna, floatrealspace> *)NULL );
}

template <typename ProfileType>
int
profillic_eslx_msafile_Read(ESLX_MSAFILE *afp, ESL_MSA **ret_msa, ProfileType * profile_ptr)
{
  ESL_MSA  *msa    = NULL;
  int       status = eslOK;
#ifdef eslAUGMENT_SSI
  esl_pos_t offset = esl_buffer_GetOffset(afp->bf);
#endif

  switch (afp->format) {
  case eslMSAFILE_A2M:          if ((status = esl_msafile_a2m_Read      (afp, &msa)) != eslOK) goto ERROR; break;
  case eslMSAFILE_AFA:          if ((status = esl_msafile_afa_Read      (afp, &msa)) != eslOK) goto ERROR; break;
  case eslMSAFILE_CLUSTAL:      if ((status = esl_msafile_clustal_Read  (afp, &msa)) != eslOK) goto ERROR; break;
  case eslMSAFILE_CLUSTALLIKE:  if ((status = esl_msafile_clustal_Read  (afp, &msa)) != eslOK) goto ERROR; break;
  case eslMSAFILE_PFAM:         if ((status = esl_msafile_stockholm_Read(afp, &msa)) != eslOK) goto ERROR; break;
  case eslMSAFILE_PHYLIP:       if ((status = esl_msafile_phylip_Read   (afp, &msa)) != eslOK) goto ERROR; break;
  case eslMSAFILE_PHYLIPS:      if ((status = esl_msafile_phylip_Read   (afp, &msa)) != eslOK) goto ERROR; break;
  case eslMSAFILE_PSIBLAST:     if ((status = esl_msafile_psiblast_Read (afp, &msa)) != eslOK) goto ERROR; break;
  case eslMSAFILE_SELEX:        if ((status = esl_msafile_selex_Read    (afp, &msa)) != eslOK) goto ERROR; break;
  case eslMSAFILE_STOCKHOLM:    if ((status = esl_msafile_stockholm_Read(afp, &msa)) != eslOK) goto ERROR; break;
  case eslMSAFILE_PROFILLIC:    if ((status = profillic_esl_msafile_profile_Read(afp, &msa, profile_ptr)) != eslOK) goto ERROR; break;
  default:                      ESL_EXCEPTION(eslEINCONCEIVABLE, "no such msa file format"); break;
  }
  
#ifdef eslAUGMENT_SSI
  msa->offset = offset;
#endif
  *ret_msa = msa;
  return eslOK;

 ERROR:
  if (msa) esl_msa_Destroy(msa);
  *ret_msa = NULL;
  return status;
}

/*****************************************************************
 * 12.5. galosh profile format (from profilic)
 *****************************************************************/

/**
 * <pre>
 *
 * Function:  profillic_esl_msafile_profile_Read()
 *
 * Paul T Edlefsen   paul@galosh.org   February 19, 2012.
 *
 * Synopsis:  Read a profillic/galosh profile.
 *
 * Purpose: Parse the Profile HMM from an open galosh profile format
 *            file (from profillic) <afp>, leaving the profile in
 *            <ret_profile>. Also create a new
 *            MSA, and return it by reference through 
 *            <*ret_msa>. Caller is responsible for freeing
 *            this <ESL_MSA>.
 *
 * Args:      <afp>     - open <ESL_MSAFILE> to read from
 *            <ret_msa> - RETURN: newly parsed, created <ESL_MSA>
 *
 * Returns:   <eslOK> on success. <*ret_msa> contains the newly
 *            allocated MSA. <afp> is poised at start of next
 *            alignment record, or is at EOF. The profile is in
 *            <ret_profile>.
 *
 *            <eslEOF> if no (more) profile data are found in
 *            <afp>, and <afp> is returned at EOF. 
 *
 *            <eslEFORMAT> on a parse error. <*ret_msa> is set to
 *            <NULL>, and <ret_profile> is unaffected.
 *            <afp> contains information sufficient for
 *            constructing useful diagnostic output: 
 *            | <afp->errmsg>       | user-directed error message     |
 *            | <afp->bf->filename> | name of the file                |
 *
 * Throws:    <eslEMEM> on allocation error.
 *            <eslESYS> if a system call fails, such as fread().
 *            <*ret_msa> is returned <NULL>.
 *</pre>      
 *      
 */
template <typename ProfileType>
static int
profillic_esl_msafile_profile_Read(ESLX_MSAFILE *afp, ESL_MSA **ret_msa, ProfileType * profile_ptr )
{
  /// \note Right now this isn't actually using the open file pointer; for convenience I just use the profile.fromFile( <filename> ) method.
  /// \todo Use convenience fns in esl_buffer.h; see eg hmmer-3.1/easel/esl_msafile_stockholm.c for examples...
  ESL_MSA                 *msa      = NULL;
  string profile_string;
  char *buf;
  long len;
  int                      seqidx;
  int                      status;
  char       errmsg2[eslERRBUFSIZE];

  ESL_DASSERT1((afp->format == eslMSAFILE_PROFILLIC));

  const char * const seqname = "Galosh Profile Consensus";
  const char * const msaname = "Galosh Profile";
  uint32_t profile_length;
  galosh::Sequence<typename ProfileType::ProfileResidueType> consensus_sequence;
  stringstream tmp_consensus_output_stream;

  uint32_t pos_i;

  if (profile_ptr == NULL)  { ESL_EXCEPTION(eslEINCONCEIVABLE, "profile_ptr is NULL in profillic_esl_msafile_profile_Read(..)!"); }
  //if (feof(afp->bf->fp))  { status = eslEOF; goto ERROR; }
  afp->errmsg[0] = '\0';

  // Read in the galosh profile (from profillic)
  //fseek( afp->bf->fp, 0, SEEK_END ); // go to the end
  //len = afp->bf->ftell( afp->bf->fp ); // get the position at the end (length)
  //fseek( afp->bf->fp, 0, SEEK_SET ); // go to the beginning again.

  //ESL_ALLOC_CPP( char, buf, sizeof( char ) * len ); //malloc buffer
  //fread( buf, len, 1, afp->bf->fp ); //read into buffer

  //profile_string = buf;
  //profile_ptr->fromString( profile_string );
  profile_ptr->fromFile( afp->bf->filename );
  //if (buf)      free(buf);
  // \todo WHY WON'T THIS WORK?  See HACKs in profillic-hmmbuild.cpp to work around it.
  //fseek( afp->bf->fp, 0, SEEK_END ); // go to the end (to signal there's no more profiles in the file, the next time we come to this function)

  // Calculate the consensus sequence.
  profile_length = profile_ptr->length();
  consensus_sequence.reinitialize( profile_length );
  for( pos_i = 0; pos_i < profile_length; pos_i++ ) {
    consensus_sequence[ pos_i ] =
      ( *profile_ptr )[ pos_i ][ galosh::Emission::Match ].maximumValueType();
  }
  tmp_consensus_output_stream << consensus_sequence;

  /* Allocate a growable MSA, and auxiliary parse data coupled to the MSA allocation */
#ifdef eslAUGMENT_ALPHABET
  if (afp->abc   &&  (msa = esl_msa_CreateDigital(afp->abc, 16, -1)) == NULL) { status = eslEMEM; goto ERROR; }
#endif
  if (! afp->abc &&  (msa = esl_msa_Create(                 16, -1)) == NULL) { status = eslEMEM; goto ERROR; }


  // Set first-and-only seq to the consensus.  This should set sqlen[0] to the profile's length and set ax to have length 1 and ax[0] to be the sequence itself.  Also msa->sqname[0] to the "name" of that consensus sequence.

  /* if nec, make room for the new seq */
  if (msa->nseq >= msa->sqalloc && (status = esl_msa_Expand(msa)) != eslOK) return status; 
  seqidx = msa->nseq; // 0
  msa->nseq++; // = 1
  status = esl_strdup(seqname, -1, &(msa->sqname[seqidx]));
  // NOTE: Could add description of this "sequence" here, using esl_msa_SetSeqDescription(msa, seqidx, desc).
#ifdef eslAUGMENT_ALPHABET
  if (msa->flags & eslMSA_DIGITAL)
    {
      // NOTE (profillic): There was a bug in this; it had said .."esl_abc_dsqcat(msa->abc, " where it should have said .."esl_abc_dsqcat(msa->abc->inmap, "
      if((status = esl_abc_dsqcat(msa->abc->inmap, &(msa->ax[seqidx]), &(msa->sqlen[seqidx]), tmp_consensus_output_stream.str().c_str(), profile_length)) != eslOK) {
        /* invalid char(s), get informative error message */
        if (esl_abc_ValidateSeq(msa->abc, tmp_consensus_output_stream.str().c_str(), profile_length, afp->errmsg) != eslOK) 
          ESL_XFAIL(eslEFORMAT, errmsg2, "%s (line %d): %s", msa->sqname[0], afp->linenumber, afp->errmsg);
      }
    }
#endif
  if (! (msa->flags & eslMSA_DIGITAL))
    {
      status = esl_strcat(&(msa->aseq[seqidx]), 0, tmp_consensus_output_stream.str().c_str(), profile_length);
      msa->sqlen[seqidx] = profile_length;
    } 
  msa->alen = profile_length;

  /// \todo OR read in a fasta file of sequences too.
  /// \todo (Optional?) Set msa->name to the name of the profile (file?)
  esl_strdup(msaname, -1, &(msa->name));
  /// \todo make sure eslMSA_HASWGTS is FALSE .. OR set it to TRUE and set msa->wgt[idx] to 1.0.
  /// \note Could have secondary structure (per sequence) too. msa->ss[0]. msa->sslen[0] should be the same as msa->sqlen[0].
  /// \todo Investigate what msa->sa and msa->pp are for.

  /* Give the newly parsed MSA a good
   * going-over, and finalize the fields of the MSA data structure.
   * verify_parse will fill in errmsg if it sees a problem.
   */
  //if (verify_parse(msa, afp->errmsg) != eslOK) { status = eslEFORMAT; goto ERROR; } 

  if (( status = esl_msa_SetDefaultWeights(msa)) != eslOK) goto ERROR;

  if (ret_msa != NULL) *ret_msa = msa; else esl_msa_Destroy(msa);
  return eslOK;

 ERROR:
  if (msa != NULL)      esl_msa_Destroy(msa);
  if (ret_msa != NULL) *ret_msa = NULL;
  return status;
}

/*---------------------- end, galosh profile format (from profillic)-------*/

/**
 * \par Licence:
 *****************************************************************
 * @LICENSE@
 *
 * SVN $URL$
 * SVN $Id: $
 *****************************************************************/

#endif // __GALOSH_PROFILLICESLMSA_HPP__
