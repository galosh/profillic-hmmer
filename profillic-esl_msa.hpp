#ifndef __GALOSH_PROFILLICESLMSA_HPP__
#define __GALOSH_PROFILLICESLMSA_HPP__

/*::cexcerpt::header_example::begin::*/
/* Multiple sequence alignment file i/o.
 *    
 * Contents:   
 *    1. The <ESL_MSA> object
 *    2. The <ESL_MSAFILE> object
 *    3. Digital mode MSA's         (augmentation: alphabet)
 *    4. Random MSA database access (augmentation: ssi)
 *    5. General i/o API, for all alignment formats
 *    6. Miscellaneous functions for manipulating MSAs
 *    7. Stockholm (Pfam/Rfam) format
 *    8. A2M format
 *    9. PSIBLAST format
 *   10. SELEX format
 *   11. AFA (aligned FASTA) format
 *   12. Memory efficient routines for PFAM format
 *   13. Debugging/development routines
 *   14. Benchmark driver
 *   15. Unit tests
 *   16. Test driver
 *   17. Examples
 *   18. Copyright and license information
 *   
 * Augmentations:
 *   alphabet:  adds support for digital MSAs
 *   keyhash:   speeds up Stockholm file input
 *   ssi:       enables indexed random access in a file of many MSAs
 *
 * to do: SRE, Sat Jan  3 09:43:42 2009 (after selex parser added)
 * - SELEX parser is better in some respects than older Stockholm
 *    parser; stricter, better error detection, better modularity.  
 *    Generalize the SELEX parser routines and use them for Stockholm.
 * - Test files for SELEX parser are in esl_msa_testfiles/selex, with
 *    tabular summary list in 00MANIFEST. This is an experiment with
 *    writing tests that require lots of external files, such as
 *    format parsers. Write test driver routine that reads 00MANIFEST
 *    and runs esl_msa_Read() against these files, checking for proper
 *    return status, including errors.
 * - The selex parser's read_block() reads lines into memory and
 *    parses them later. afp->linenumber is thus no longer an
 *    accurate record of where a parse error occurs. read_xxx()
 *    format parsers now need to include line number in their 
 *    afp->errbuf[] message upon eslEFORMAT error. Stockholm parser
 *    doesn't do this. Make it so, and document in examples.
 * - Format autodetection doesn't work yet. Coordinate w/ how sqio
 *    does it, and implement. May require buffering input to make
 *    it work with .gz, pipes without rewinding a stream. Might be
 *    a good idea to generalize input buffering - perhaps making
 *    it part of ESL_FILEPARSER. 
 * - PSIBLAST, A2M format only supported on output, not input.
 *    Implement input parsers.
 * - SELEX format only supported on input, not output. 
 *    Implement output writer.
 * - More formats need to be parsed. Check on formats for current
 *    best MSA programs, such as MUSCLE, MAFFT; implement i/o.
 *    
 * SRE, Thu Jan 20 08:50:43 2005 [St. Louis]
 * SVN $Id: esl_msa.c 573 2010-03-27 15:13:52Z eddys $
 */
/*::cexcerpt::header_example::end::*/

/*::cexcerpt::include_example::begin::*/
extern "C" {
#include "esl_config.h"
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

extern "C" {
#include "easel.h"
#ifdef eslAUGMENT_KEYHASH
#include "esl_keyhash.h"
#endif
#ifdef eslAUGMENT_ALPHABET
#include "esl_alphabet.h"
#endif
#ifdef eslAUGMENT_SSI
#include "esl_ssi.h"
#endif
#include "esl_msa.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"
}
/*::cexcerpt::include_example::end::*/



/////////////// For profillic-hmmer //////////////////////////////////
#include "profillic-hmmer.hpp"
extern "C" {
#include "esl_msa.h"
}
#define eslMSAFILE_PROFILLIC       98103  /* A galosh profile (from profillic)   */
/////////////// End profillic-hmmer //////////////////////////////////

/******************************************************************************
 *# 1. The <ESL_MSA> object                                           
 *****************************************************************************/
/* get_seqidx()
 * 
 * Find the index of a given sequence <name> in an <msa>.
 * If caller has a good guess (for instance, the sequences are
 * coming in a previously seen order in a block of seqs or annotation),
 * the caller can pass this information in <guess>, or -1 if
 * it has no guess.
 * 
 * This function behaves differently depending on whether
 * keyhash augmentation is available or not. Without keyhashing,
 * the name is identified by bruteforce search of the names
 * in the <msa>. With keyhashing, we hash search, which should
 * improve performance for large alignments.
 * 
 * If the name does not already exist in the MSA, then it
 * is assumed to be a new sequence name that we need to store.
 * seqidx is set to msa->nseq, the MSA is Expand()'ed if necessary
 * to make room, the name is stored in msa->sqname[msa->nseq],
 * (and in the hash table, if we're keyhash augmented)
 * and msa->nseq is incremented.
 *
 * Returns:  <eslOK> on success, and the seqidx is 
 *           passed back via <ret_idx>. If <name> is new
 *           in the <msa>, the <name> is stored and the <msa> 
 *           may be internally reallocated if needed.
 *           
 * Throws:   <eslEMEM> if we try to add a name and allocation fails.
 *           <eslEINVAL> if we try to add a name to a non-growable MSA.
 */
static int
get_seqidx(ESL_MSA *msa, char *name, int guess, int *ret_idx)
{
  int seqidx;
  int status;

  *ret_idx = -1;

  /* can we guess? */
  if (guess >= 0 && 
      guess < msa->nseq && 
      strcmp(name, msa->sqname[guess]) == 0) 
    { *ret_idx = guess; return eslOK; }

  /* Else look it up - either brute force
   * or, if we're keyhash-augmented, by hashing.
   */
#ifdef eslAUGMENT_KEYHASH                  
  status = esl_key_Store(msa->index, name, &seqidx);
  if (status == eslEDUP) { *ret_idx = seqidx; return eslOK; }
  if (status != eslOK) return status; /* an error. */
#else
  for (seqidx = 0; seqidx < msa->nseq; seqidx++)
    if (strcmp(msa->sqname[seqidx], name) == 0) break;
  if (seqidx < msa->nseq) 
    { *ret_idx = seqidx; return eslOK; }
#endif

  /* If we reach here, then this is a new name that we're
   * adding.
   */
  if (seqidx >= msa->sqalloc &&  
     (status = esl_msa_Expand(msa)) != eslOK)
    return status; 
    
  status = esl_strdup(name, -1, &(msa->sqname[seqidx]));
  msa->nseq++;
  if (ret_idx != NULL) *ret_idx = seqidx;
  return status;
}



/* verify_parse()
 *
 * Last function called after a multiple alignment parser thinks it's
 * done. Checks that parse was successful; makes sure required
 * information is present; makes sure required information is
 * consistent. Some fields that are only use during parsing may be
 * freed (sqlen, for example), and some fields are finalized now
 * (<msa->alen> is set, for example). 
 * 
 * <errbuf> is a place to sprintf an informative message about the
 * reason for a parse error. The caller provides an <errbuf>
 * of at least 512 bytes.
 *
 * Returns:  <eslOK>, and errbuf is set to an empty string.
 *           
 * Throws:   <eslEFORMAT> if a problem is detected, and an
 *           informative message about the failure is in errbuf.
 */
static int
verify_parse(ESL_MSA *msa, char *errbuf)
{
  int idx;

  if (msa->nseq == 0) ESL_FAIL(eslEFORMAT, errbuf, "parse error: no alignment data found");

  /* set alen, until proven otherwise; we'll check that the other seqs
   * have the same length later.
   */
  msa->alen = msa->sqlen[0];

  /* We can rely on msa->sqname[] being valid for any index,
   * because of the way the line parsers always store any name
   * they add to the index.
   */
  for (idx = 0; idx < msa->nseq; idx++)
    {
#ifdef eslAUGMENT_ALPHABET
      if ((msa->flags & eslMSA_DIGITAL) &&  (msa->ax  == NULL || msa->ax[idx] == NULL))
	ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: no sequence for %s",
		 msa->name != NULL ? msa->name : "", msa->sqname[idx]); 
#endif
      if (! (msa->flags & eslMSA_DIGITAL) && (msa->aseq == NULL || msa->aseq[idx] == NULL))
	ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: no sequence for %s",
		 msa->name != NULL ? msa->name : "", msa->sqname[idx]); 

      /* either all weights must be set, or none of them */
      if ((msa->flags & eslMSA_HASWGTS) && msa->wgt[idx] == -1.0)
	ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: expected a weight for seq %s", 
		  msa->name != NULL ? msa->name : "", msa->sqname[idx]);

      /* all aseq must be same length. */
      if (msa->sqlen[idx] != msa->alen)
	ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: sequence %s: length %ld, expected %ld",
		 msa->name != NULL ? msa->name : "", msa->sqname[idx], static_cast<long int>( msa->sqlen[idx] ), static_cast<long int>( msa->alen ) );

      /* if individual SS is present, it must have length right too */
      if (msa->ss != NULL &&  msa->ss[idx] != NULL &&  msa->sslen[idx] != msa->alen) 
	ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: GR SS for %s: length %ld, expected %ld",
		 msa->name != NULL ? msa->name : "", msa->sqname[idx], msa->sslen[idx], msa->alen);

				/* if SA is present, must have length right */
      if (msa->sa != NULL && msa->sa[idx] != NULL && msa->salen[idx] != msa->alen) 
	ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: GR SA for %s: length %ld, expected %ld",
		 msa->name != NULL ? msa->name : "", msa->sqname[idx], msa->salen[idx], msa->alen);

				/* if PP is present, must have length right */
      if (msa->pp != NULL && msa->pp[idx] != NULL && msa->pplen[idx] != msa->alen) 
	ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: GR PP for %s: length %ld, expected %ld",
		 msa->name != NULL ? msa->name : "", msa->sqname[idx], msa->pplen[idx], msa->alen);
    }

  /* if cons SS is present, must have length right */
  if (msa->ss_cons != NULL && strlen(msa->ss_cons) != msa->alen) 
    ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: GC SS_cons markup: len %zd, expected %ld",
	     msa->name != NULL ? msa->name : "",  strlen(msa->ss_cons), msa->alen);

  /* if cons SA is present, must have length right */
  if (msa->sa_cons != NULL && strlen(msa->sa_cons) != msa->alen) 
    ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: GC SA_cons markup: len %zd, expected %ld",
	     msa->name != NULL ? msa->name : "",  strlen(msa->sa_cons), msa->alen);

  /* if cons PP is present, must have length right */
  if (msa->pp_cons != NULL && strlen(msa->pp_cons) != msa->alen) 
    ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: GC PP_cons markup: len %zd, expected %ld",
	     msa->name != NULL ? msa->name : "",  strlen(msa->pp_cons), msa->alen);

  /* if RF is present, must have length right */
  if (msa->rf != NULL && strlen(msa->rf) != msa->alen) 
    ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: GC RF markup: len %zd, expected %ld",
	     msa->name != NULL ? msa->name : "", strlen(msa->rf), msa->alen);

  /* If no weights were set, set 'em all to 1.0 */
  if (!(msa->flags & eslMSA_HASWGTS))
    for (idx = 0; idx < msa->nseq; idx++)
      msa->wgt[idx] = 1.0;

  /* Clean up a little from the parser */
  if (msa->sqlen != NULL) { free(msa->sqlen); msa->sqlen = NULL; }
  if (msa->sslen != NULL) { free(msa->sslen); msa->sslen = NULL; }
  if (msa->salen != NULL) { free(msa->salen); msa->salen = NULL; }
  if (msa->pplen != NULL) { free(msa->pplen); msa->pplen = NULL; }
  return eslOK;
}

static int read_stockholm(ESL_MSAFILE *afp, ESL_MSA **ret_msa);
static int read_selex    (ESL_MSAFILE *afp, ESL_MSA **ret_msa);
static int read_afa      (ESL_MSAFILE *afp, ESL_MSA **ret_msa);
template <typename ProfileType>
static int profillic_read_profile      (ESL_MSAFILE *afp, ESL_MSA **ret_msa, ProfileType * profile_ptr);

/* Function:  esl_msa_Read()
 * Synopsis:  Read next MSA from a file.
 * Incept:    SRE, Fri Jan 28 08:10:49 2005 [St. Louis]
 *
 * Purpose:   Reads the next MSA from an open MSA file <afp>,
 *            and returns it via <ret_msa>. 
 *
 * Returns:   <eslOK> on success, and <ret_msa> points at the
 *            new MSA object.
 *            
 *            Returns <eslEOF> if there are no more alignments in the file.
 *            
 *            Returns <eslEFORMAT> if there is a parse error, and <afp->errbuf>
 *            is set to an informative message.
 *            
 *            <eslEINVAL> if we're trying to read a digital alignment,
 *            but one or more residues are seen in the file that
 *            aren't valid in our alphabet.
 *            
 * Throws:    <eslEMEM> on allocation failure.           
 *            <eslEINCONCEIVABLE> on internal error.
 */
template <class ProfileType>
int
profillic_esl_msa_Read(ESL_MSAFILE *afp, ESL_MSA **ret_msa, ProfileType * profile_ptr)
{
  ESL_MSA *msa;
  int      status;

  *ret_msa = NULL;

  /* If we've just used GuessAlphabet(), we have an MSA already read
   * and stored in the MSAFILE's cache. Just return it, after worrying
   * about whether it's supposed to be in digital or text mode. (It
   * should always be in text mode, and maybe in need of Digitize(),
   * given how GuessAlphabet works now; but this is coded for more
   * generality in case we use the MSA cache some other way in the
   * future.)
   */
  if (afp->msa_cache != NULL) 
    {
#ifdef eslAUGMENT_ALPHABET
      if      (afp->do_digital   && !(afp->msa_cache->flags & eslMSA_DIGITAL)) {
	if ((status = esl_msa_Digitize(afp->abc, afp->msa_cache, afp->errbuf)) != eslOK) return status; 
      }
      else if (! afp->do_digital && (afp->msa_cache->flags & eslMSA_DIGITAL)) {
	if ((status = esl_msa_Textize(afp->msa_cache)) != eslOK) return status;
      }
#endif

      *ret_msa         = afp->msa_cache;
      afp->msa_cache = NULL;
      return eslOK;
    }

  /* Otherwise, read the next MSA from the file.
   */      
  switch (afp->format) {
  case eslMSAFILE_STOCKHOLM: status = read_stockholm(afp, &msa); break;
  case eslMSAFILE_PFAM:      status = read_stockholm(afp, &msa); break;
  case eslMSAFILE_A2M:       ESL_FAIL(eslEFORMAT, afp->errbuf, "A2M format input parser not implemented yet.");
  case eslMSAFILE_PSIBLAST:  ESL_FAIL(eslEFORMAT, afp->errbuf, "PSIBLAST format input parser not implemented yet.");
  case eslMSAFILE_SELEX:     status = read_selex    (afp, &msa); break;
  case eslMSAFILE_AFA:       status = read_afa      (afp, &msa); break;
  case eslMSAFILE_PROFILLIC:  status = profillic_read_profile      (afp, &msa, profile_ptr); break;
  default:                   ESL_EXCEPTION(eslEINCONCEIVABLE, "no such format");
  }

  *ret_msa = msa;
  return status;
}


/* Function:  esl_msa_EncodeFormat()
 * Synopsis:  Convert text string to an MSA file format code.
 * Incept:    SRE, Fri Oct 24 13:21:08 2008 [Janelia]
 *
 * Purpose:   Given a text string, match it case-insensitively
 *            against a list of possible formats, and return the
 *            appropriate MSA file format code. For example,
 *            <esl_msa_EncodeFormat("Stockholm")> returns
 *            <eslMSAFILE_STOCKHOLM>.
 *            
 *            If the format is unrecognized, return
 *            <eslMSAFILE_UNKNOWN>.
 *            
 * Note:      Keep in sync with <esl_sqio_EncodeFormat()>, 
 *            which decodes all possible sequence file formats,
 *            both unaligned and aligned.           
 */
int
profillic_esl_msa_EncodeFormat(char *fmtstring)
{
  if (strcasecmp(fmtstring, "stockholm") == 0) return eslMSAFILE_STOCKHOLM;
  if (strcasecmp(fmtstring, "pfam")      == 0) return eslMSAFILE_PFAM;
  if (strcasecmp(fmtstring, "a2m")       == 0) return eslMSAFILE_A2M;
  if (strcasecmp(fmtstring, "psiblast")  == 0) return eslMSAFILE_PSIBLAST;
  if (strcasecmp(fmtstring, "selex")     == 0) return eslMSAFILE_SELEX;
  if (strcasecmp(fmtstring, "afa")       == 0) return eslMSAFILE_AFA;
  if (strcasecmp(fmtstring, "profillic") == 0) return eslMSAFILE_PROFILLIC;
  return eslMSAFILE_UNKNOWN;
}


/*****************************************************************
 * 7. Stockholm (Pfam/Rfam) format
 *****************************************************************/

/* msafile_getline():
 * load the next line of <afp> into <afp->buf>. 
 * Returns eslOK on success, eslEOF on normal eof.
 * Throws eslEMEM on alloc failure.
 */
static int
msafile_getline(ESL_MSAFILE *afp)
{
  int status;
  status = esl_fgets(&(afp->buf), &(afp->buflen), afp->f);
  afp->linenumber++;
  return status;
}

/* maxwidth()
 * Return the length of the longest string in 
 * an array of strings.
 */
static int64_t
maxwidth(char **s, int n)
{
  int64_t max,len;
  int     i; 
  
  max = 0;
  for (i = 0; i < n; i++)
    if (s[i] != NULL)
      {
	len = strlen(s[i]);
	if (len > max) max = len;
      }
  return max;
}
static int
is_blankline(char *s)
{
  for (; *s != '\0'; s++)
    if (! isspace((int) *s)) return FALSE;
  return TRUE;
}

/* Format of a GF line:
 *    #=GF <tag> <text>
 * Returns eslOK on success; eslEFORMAT on parse failure.
 * Throws eslEMEM on allocation failure.
 */
static int
parse_gf(ESL_MSA *msa, char *buf)
{
  char *gf;
  char *tag;
  char *text;
  char *tok;
  char *s;
  int   n;
  int   status;

  s = buf;
  if (esl_strtok(&s, " \t\n\r", &gf)  != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, " \t\n\r", &tag) != eslOK) return eslEFORMAT;

  /* text might be empty; watch out for this. (for example, a blank #=GF CC line) */
  status = esl_strtok_adv(&s, "\n\r",    &text, &n, NULL);
  if      (status == eslOK) { while (*text && (*text == ' ' || *text == '\t')) text++; }
  else if (status == eslEOL){ text = NULL; n = 0; } 
  else return eslEFORMAT;

  if      (strcmp(tag, "ID") == 0) status = esl_strdup(text, n, &(msa->name));
  else if (strcmp(tag, "AC") == 0) status = esl_strdup(text, n, &(msa->acc));
  else if (strcmp(tag, "DE") == 0) status = esl_strdup(text, n, &(msa->desc));
  else if (strcmp(tag, "AU") == 0) status = esl_strdup(text, n, &(msa->au));
  else if (strcmp(tag, "GA") == 0) 
    {				/* Pfam has GA1, GA2. Rfam just has GA1. */
      s = text;
      if ((esl_strtok(&s, " \t\n\r", &tok)) != eslOK) 
	return eslEFORMAT;
      msa->cutoff[eslMSA_GA1] = atof(tok);
      msa->cutset[eslMSA_GA1] = TRUE;
      if ((esl_strtok(&s, " \t\n\r", &tok)) == eslOK) 
	{
	  msa->cutoff[eslMSA_GA2] = atof(tok);
	  msa->cutset[eslMSA_GA2] = TRUE;
	}
      status = eslOK;
    }
  else if (strcmp(tag, "NC") == 0) 
    {
      s = text;
      if ((esl_strtok(&s, " \t\n\r", &tok)) != eslOK) 
	return eslEFORMAT;
      msa->cutoff[eslMSA_NC1] = atof(tok);
      msa->cutset[eslMSA_NC1] = TRUE;
      if ((esl_strtok(&s, " \t\n\r", &tok)) == eslOK) 
	{
	  msa->cutoff[eslMSA_NC2] = atof(tok);
	  msa->cutset[eslMSA_NC2] = TRUE;
	}
      status = eslOK;
    }
  else if (strcmp(tag, "TC") == 0) 
    {
      s = text;
      if ((esl_strtok(&s, " \t\n\r", &tok)) != eslOK) 
	return eslEFORMAT;
      msa->cutoff[eslMSA_TC1] = atof(tok);
      msa->cutset[eslMSA_TC1] = TRUE;
      if ((esl_strtok(&s, "\t\n\r", &tok)) == eslOK) 
	{
	  msa->cutoff[eslMSA_TC2] = atof(tok);
	  msa->cutset[eslMSA_TC2] = TRUE;
	}
      status = eslOK;
    }
  else 				/* an unparsed #=GF: */
    status = esl_msa_AddGF(msa, tag, text);

  return status;
}


/* Format of a GS line:
 *    #=GS <seqname> <tag> <text>
 * Return <eslOK> on success; <eslEFORMAT> on parse error.
 * Throws <eslEMEM> on allocation error (trying to grow for a new
 *        name; <eslEINVAL> if we try to grow an ungrowable MSA.
 */
static int
parse_gs(ESL_MSA *msa, char *buf)
{
  char *gs;
  char *seqname;
  char *tag;
  char *text; 
  int   seqidx;
  char *s;
  int   status;

  s = buf;
  if (esl_strtok(&s, " \t\n\r", &gs)      != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, " \t\n\r", &seqname) != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, " \t\n\r", &tag)     != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, "\n\r",    &text)    != eslOK) return eslEFORMAT;
  while (*text && (*text == ' ' || *text == '\t')) text++;
  
  /* GS usually follows another GS; guess lastidx+1 */
  status = get_seqidx(msa, seqname, msa->lastidx+1, &seqidx);
  if (status != eslOK) return status;
  msa->lastidx = seqidx;

  if (strcmp(tag, "WT") == 0)
    {
      msa->wgt[seqidx] = atof(text);
      msa->flags      |= eslMSA_HASWGTS;
      status           = eslOK;
    }
  else if (strcmp(tag, "AC") == 0)
    status = esl_msa_SetSeqAccession(msa, seqidx, text);
  else if (strcmp(tag, "DE") == 0)
    status = esl_msa_SetSeqDescription(msa, seqidx, text);
  else				
    status = esl_msa_AddGS(msa, tag, seqidx, text);

  return status;
}


/* parse_gc():
 * Format of a GC line:
 *    #=GC <tag> <aligned text>
 */
static int 
parse_gc(ESL_MSA *msa, char *buf)
{
  char *gc;
  char *tag;
  char *text; 
  char *s;
  int   len;
  int   status;

  s = buf;
  if (esl_strtok    (&s, " \t\n\r", &gc)               != eslOK) return eslEFORMAT;
  if (esl_strtok    (&s, " \t\n\r", &tag)              != eslOK) return eslEFORMAT;
  if (esl_strtok_adv(&s, " \t\n\r", &text, &len, NULL) != eslOK) return eslEFORMAT;
  
  if      (strcmp(tag, "SS_cons") == 0)  status = esl_strcat(&(msa->ss_cons), -1, text, len);
  else if (strcmp(tag, "SA_cons") == 0)  status = esl_strcat(&(msa->sa_cons), -1, text, len);
  else if (strcmp(tag, "PP_cons") == 0)  status = esl_strcat(&(msa->pp_cons), -1, text, len);
  else if (strcmp(tag, "RF")      == 0)  status = esl_strcat(&(msa->rf),      -1, text, len);
  else                                   status = esl_msa_AppendGC(msa, tag, text);

  return status;
}

/* parse_gr():
 * Format of a GR line:
 *    #=GR <seqname> <featurename> <text>
 */
static int
parse_gr(ESL_MSA *msa, char *buf)
{
  char *gr;
  char *seqname;
  char *tag;
  char *text;
  int   seqidx;
  int   len;
  int   j;
  char *s;
  int   status;

  s = buf;
  if (esl_strtok    (&s, " \t\n\r", &gr)               != eslOK) return eslEFORMAT;
  if (esl_strtok    (&s, " \t\n\r", &seqname)          != eslOK) return eslEFORMAT;
  if (esl_strtok    (&s, " \t\n\r", &tag)              != eslOK) return eslEFORMAT;
  if (esl_strtok_adv(&s, " \t\n\r", &text, &len, NULL) != eslOK) return eslEFORMAT;

  /* GR usually follows sequence it refers to; guess msa->lastidx */
  status = get_seqidx(msa, seqname, msa->lastidx, &seqidx);
  if (status != eslOK) return status;
  msa->lastidx = seqidx;

  if (strcmp(tag, "SS") == 0) 
    {
      if (msa->ss == NULL)
	{
	  ESL_ALLOC_CPP(char *, msa->ss,    sizeof(char *) * msa->sqalloc);
	  ESL_ALLOC_CPP(int64_t, msa->sslen, sizeof(int64_t)* msa->sqalloc);
	  for (j = 0; j < msa->sqalloc; j++)
	    {
	      msa->ss[j]    = NULL;
	      msa->sslen[j] = 0;
	    }
	}
      status = esl_strcat(&(msa->ss[seqidx]), msa->sslen[seqidx], text, len);
      msa->sslen[seqidx] += len;
    }
  else if (strcmp(tag, "SA") == 0)
    {
      if (msa->sa == NULL)
	{
	  ESL_ALLOC_CPP(char *, msa->sa,    sizeof(char *) * msa->sqalloc);
	  ESL_ALLOC_CPP(int64_t, msa->salen, sizeof(int64_t)* msa->sqalloc);
	  for (j = 0; j < msa->sqalloc; j++) 
	    {
	      msa->sa[j]    = NULL;
	      msa->salen[j] = 0;
	    }
	}
      status = esl_strcat(&(msa->sa[seqidx]), msa->salen[seqidx], text, len);
      msa->salen[seqidx] += len;
    }
  else if (strcmp(tag, "PP") == 0)
    {
      if (msa->pp == NULL)
	{
	  ESL_ALLOC_CPP(char *, msa->pp,    sizeof(char *) * msa->sqalloc);
	  ESL_ALLOC_CPP(int64_t, msa->pplen, sizeof(int64_t)* msa->sqalloc);
	  for (j = 0; j < msa->sqalloc; j++) 
	    {
	      msa->pp[j]    = NULL;
	      msa->pplen[j] = 0;
	    }
	}
      status = esl_strcat(&(msa->pp[seqidx]), msa->pplen[seqidx], text, len);
      msa->pplen[seqidx] += len;
    }
  else 
    status = esl_msa_AppendGR(msa, tag, seqidx, text);
  return status;

 ERROR:
  return status;
}


/* parse_comment():
 * comments are simply stored verbatim, not parsed
 */
static int
parse_comment(ESL_MSA *msa, char *buf)
{
  char *s;
  char *comment;

  s = buf + 1;			               /* skip leading '#' */
  if (*s == '\n' || *s == '\r') { *s = '\0'; comment = s; }  /* deal with blank comment */
  else if (esl_strtok(&s, "\n\r", &comment)!= eslOK) return eslEFORMAT;
  return (esl_msa_AddComment(msa, comment));
}

/* parse_sequence():
 * Format of line is:
 *     <name>  <aligned text>
 * 
 * On digital sequence, returns <eslEINVAL> if any of the residues can't be digitized.
 */
static int
parse_sequence(ESL_MSA *msa, char *buf)
{
  char *s;
  char *seqname;
  char *text;
  int   seqidx;
  int   len;
  int   status;

  s = buf;
  if (esl_strtok    (&s, " \t\n\r", &seqname)          != eslOK) return eslEFORMAT;
  if (esl_strtok_adv(&s, " \t\n\r", &text, &len, NULL) != eslOK) return eslEFORMAT; 
  
  /* seq usually follows another seq; guess msa->lastidx +1 */
  status = get_seqidx(msa, seqname, msa->lastidx+1, &seqidx);
  if (status != eslOK) return status;
  msa->lastidx = seqidx;

#ifdef eslAUGMENT_ALPHABET
  if (msa->flags & eslMSA_DIGITAL)
    {
      status = esl_abc_dsqcat(msa->abc, &(msa->ax[seqidx]), &(msa->sqlen[seqidx]), text, len);
    }
#endif
  if (! (msa->flags & eslMSA_DIGITAL))
    {
      status = esl_strcat(&(msa->aseq[seqidx]), msa->sqlen[seqidx], text, len);
      msa->sqlen[seqidx] += len;
    }

  return status;
}

/* read_stockholm():
 * SRE, Sun Jan 23 08:33:32 2005 [St. Louis]
 *
 * Purpose:   Parse the next alignment from an open Stockholm format alignment
 *            file <afp>, leaving the alignment in <ret_msa>.
 *
 * Returns:   <eslOK> on success, and the alignment is in <ret_msa>.
 *            Returns <eslEOF> if there are no more alignments in <afp>,
 *            and <ret_msa> is set to NULL.
 *            <eslEFORMAT> if parse fails because of a file format problem,
 *            in which case afp->errbuf is set to contain a formatted message 
 *            that indicates the cause of the problem, and <ret_msa> is
 *            set to NULL. 
 *
 *            Returns <eslEINVAL> if we're trying to read a digital alignment,
 *            and an invalid residue is found that can't be digitized.
 *
 * Throws:    <eslEMEM> on allocation error.
 *
 * Xref:      squid's ReadStockholm(), 1999.
 */
static int
read_stockholm(ESL_MSAFILE *afp, ESL_MSA **ret_msa)
{
  ESL_MSA   *msa = NULL;
  char      *s;
  int        status;
  int        status2;
#ifdef eslAUGMENT_SSI
  off_t      offset;
#endif

  if (feof(afp->f))  { status = eslEOF; goto ERROR; }
  afp->errbuf[0] = '\0';

  /* Initialize allocation of the MSA:
   * make it growable, by giving it an initial blocksize of
   * 16 seqs of 0 length.
   */
#ifdef eslAUGMENT_ALPHABET
  if (afp->do_digital == TRUE && (msa = esl_msa_CreateDigital(afp->abc, 16, -1))  == NULL) 
    { status = eslEMEM; goto ERROR; }

#endif
  if (afp->do_digital == FALSE && (msa = esl_msa_Create(16, -1))  == NULL)
    { status = eslEMEM; goto ERROR; }
  if (msa == NULL)    
    { status = eslEMEM; goto ERROR; }

  /* Check the magic Stockholm header line.
   * We have to skip blank lines here, else we perceive
   * trailing blank lines in a file as a format error when
   * reading in multi-record mode.
   */
  do {
#ifdef eslAUGMENT_SSI
    offset = ftello(afp->f);
#endif
    if ((status = msafile_getline(afp)) != eslOK) goto ERROR; /* includes EOF  */
  } while (is_blankline(afp->buf));

  if (strncmp(afp->buf, "# STOCKHOLM 1.", 14) != 0)
    ESL_XFAIL(eslEFORMAT, afp->errbuf, "parse failed (line %d): missing \"# STOCKHOLM\" header", afp->linenumber);

#ifdef eslAUGMENT_SSI
  msa->offset = offset;
#endif

  /* Read the alignment file one line at a time.
   */
  while ((status2 = msafile_getline(afp)) == eslOK) 
    {
      s = afp->buf;
      while (*s == ' ' || *s == '\t') s++;  /* skip leading whitespace */

      if (*s == '#') {

	if      (strncmp(s, "#=GF", 4) == 0)
	  {
	    if ((status = parse_gf(msa, s)) != eslOK)
	      ESL_XFAIL(status, afp->errbuf, "parse failed (line %d): bad #=GF line", afp->linenumber);
	  }

	else if (strncmp(s, "#=GS", 4) == 0)
	  {
	    if ((status = parse_gs(msa, s)) != eslOK)
	      ESL_XFAIL(status, afp->errbuf, "parse failed (line %d): bad #=GS line", afp->linenumber);
	  }

	else if (strncmp(s, "#=GC", 4) == 0)
	  {
	    if  ((status = parse_gc(msa, s)) != eslOK)
	      ESL_XFAIL(status, afp->errbuf, "parse failed (line %d): bad #=GC line", afp->linenumber);
	  }

	else if (strncmp(s, "#=GR", 4) == 0)
	  {
	    if ((status = parse_gr(msa, s)) != eslOK)
	      ESL_XFAIL(status, afp->errbuf, "parse failed (line %d): bad #=GR line", afp->linenumber);
	  }

	else if ((status = parse_comment(msa, s)) != eslOK)
	  ESL_XFAIL(status, afp->errbuf, "parse failed (line %d): bad comment line", afp->linenumber);
      } 
      else if (strncmp(s, "//",   2) == 0)   break; /* normal way out */
      else if (*s == '\n' || *s == '\r')     continue;
      else if ((status = parse_sequence(msa, s)) != eslOK)
	ESL_XFAIL(status, afp->errbuf, "parse failed (line %d): bad sequence line", afp->linenumber);
    }
  /* If we saw a normal // end, we would've successfully read a line,
   * so when we get here, status (from the line read) should be eslOK.
   */ 
  if (status2 != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "parse failed (line %d): didn't find // at end of alignment", afp->linenumber);
  
  /* Stockholm fmt is complex, so give the newly parsed MSA a good
   * going-over, and finalize the fields of the MSA data structure.
   * verify_parse will fill in errbuf if it sees a problem.
   */
  if (verify_parse(msa, afp->errbuf) != eslOK) { status = eslEFORMAT; goto ERROR; } 

  if (ret_msa != NULL) *ret_msa = msa; else esl_msa_Destroy(msa);
  return eslOK;

 ERROR:
  if (msa != NULL)      esl_msa_Destroy(msa);
  if (ret_msa != NULL) *ret_msa = NULL;
  return status;
}


/*****************************************************************
 * 10. SELEX format
 *****************************************************************/
#define eslMSA_LINE_SQ 1
#define eslMSA_LINE_RF 2
#define eslMSA_LINE_CS 3
#define eslMSA_LINE_SS 4
#define eslMSA_LINE_SA 5

static int read_block(ESL_MSAFILE *afp, char ***line_p, int **llen_p, int **lpos_p, int **rpos_p, int *lalloc_p, int *nlines_p, int *ret_starti);
static int first_selex_block(char *errbuf, int starti, char **line, int *lpos, int *rpos, int nlines, ESL_MSA **ret_msa, int **ret_ltype);
static int other_selex_block(char *errbuf, int starti, char **line, int *lpos, int *rpos, int nlines, ESL_MSA      *msa, int      *ltype);
static int append_selex_block(ESL_MSA *msa, char **line, int *ltype, int *lpos, int *rpos, int nlines);

/* read_selex()
 * Read an alignment in SELEX format.
 * SRE, Mon Dec 29 10:19:32 2008 [Pamplona]
 *
 * Purpose:  Parse an alignment from an open SELEX format alignment
 *           file <afp>, returning the alignment in <ret_msa>.
 *
 * Returns:   <eslOK> on success, and the alignment is in <ret_msa>.
 *
 *            Returns <eslEFORMAT> if parse fails because of a file
 *            format problem. 
 *            Returns <eslEOF> if no alignment is found in the file.
 *            Returns <eslEINVAL> if we're trying to read a digital
 *            alignment, and an invalid residue is found that 
 *            can't be digitized.
 *
 *            On all normal error conditions, <afp->errbuf> contains
 *            an informative error message for the user, and the 
 *            <*ret_msa> is <NULL>. The error message looks like
 *            "parse failed (line 156): too many #=SS lines for seq"
 *            The caller can prefix with filename if it likes.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
static int
read_selex(ESL_MSAFILE *afp, ESL_MSA **ret_msa)
{
  ESL_MSA  *msa     = NULL;
  char    **line    = NULL;
  int      *ltype   = NULL;
  int      *llen    = NULL;
  int      *lpos    = NULL;
  int      *rpos    = NULL;
  int       lalloc  = 0;
  int       nlines  = 0;
  int       nblocks = 0;
  int       starti;
  int       i, apos;
  int       status;

  if (feof(afp->f))  { status = eslEOF; goto ERROR; }
  afp->errbuf[0] = '\0';

  /* For each alignment block: */
  while ( (status = read_block(afp, &line, &llen, &lpos, &rpos, &lalloc, &nlines, &starti)) == eslOK)
    { /* now line[0..nlines-1] are data lines; llen[0..nlines-1] are max line lengths exc. \0 */
      /* lpos[0..nlines-1] are 0; and rpos[0..nlines-1] are idx of last nonwhitespace char on lines */
      nblocks++;

      if (nblocks == 1) status = first_selex_block(afp->errbuf, starti, line, lpos, rpos, nlines, &msa, &ltype);
      else              status = other_selex_block(afp->errbuf, starti, line, lpos, rpos, nlines,  msa,  ltype);
      if (status != eslOK) goto ERROR;

      if ((status = append_selex_block(msa, line, ltype, lpos, rpos, nlines)) != eslOK) goto ERROR;
    }
  if (status != eslEOF || nblocks == 0) goto ERROR; 

#ifdef eslAUGMENT_SSI
  msa->offset = 0; /* SELEX files are single MSA only; offset is always 0. */
#endif

  /* SELEX format allows ' ' as gaps, but easel doesn't */
  if (msa->rf != NULL)  
    for (apos = 0; apos < msa->alen; apos++)
      if (msa->rf[apos] == ' ') msa->rf[apos] = '.';
  if (msa->ss_cons != NULL)  
    for (apos = 0; apos < msa->alen; apos++)
      if (msa->ss_cons[apos] == ' ') msa->ss_cons[apos] = '.';
  if (msa->ss != NULL) 
    for (i = 0; i < msa->nseq; i++)
      if (msa->ss[i] != NULL) 
	for (apos = 0; apos < msa->alen; apos++)
	  if (msa->ss[i][apos] == ' ') msa->ss[i][apos] = '.';
  if (msa->sa != NULL) 
    for (i = 0; i < msa->nseq; i++)
      if (msa->sa[i] != NULL) 
	for (apos = 0; apos < msa->alen; apos++)
	  if (msa->sa[i][apos] == ' ') msa->sa[i][apos] = '.';
  for (i = 0; i < msa->nseq; i++)
    for (apos = 0; apos < msa->alen; apos++)
      if (msa->aseq[i][apos] == ' ') msa->aseq[i][apos] = '.';

#ifdef eslAUGMENT_ALPHABET 
  if (afp->do_digital) status = esl_msa_Digitize(afp->abc, msa, afp->errbuf);
#endif  
    
  *ret_msa = msa;
  free(ltype);
  return eslOK;
  
 ERROR:
  if (msa     != NULL) esl_msa_Destroy(msa);
  if (ltype   != NULL) free(ltype); 
  return status;
}

/* read_block()
 * Read one block of alignment data into memory.
 *
 * If we're trying to read the *first* block in an alignment, then at start:
 *   No lines of the alignment file have been read into <afp->buf> yet. 
 *   *line_p, *llen_p, *lpos_p, *rpos_p are all NULL
 *   *lalloc_p is 0
 *   *nlines_p is 0
 * On success, returns <eslOK> (even if last line is end of file), and:
 *   afp->buf either contains a blank line (immediately after block end), or <afp> is at EOF.
 *   *nlines_p points to the number of lines stored
 *   *line_p points to line[0..nlines-1][] array of \0 terminated strings
 *   *llen_p points to llen[0..nlines-1] array of string allocations in chars, not including \0
 *   *lpos_p points to lpos[0..nlines-1] array, all initialized to 0
 *   *rpos_p points to rpos[0..nlines-1] array, all initialized to idx of last nonwhitespace char on line
 *   *lalloc_p is >= *nlines
 * If file is empty (no data), returns <eslEOF>.
 * If an allocation fails at any point, throw <eslEMEM>. 
 *
 * If we are trying to read a *subsequent* block in the file, at start:
 *   <afp->buf> is where a previous read left it: on a blank line, or <afp> is at EOF.
 *   *line_p points to previous line[][] array; we'll reuse it, reallocating if needed.
 *   *llen_p points to previous llen[] array; ditto
 *   *lpos_p points to lpos[] array: all 0, same as first block
 *   *rpos_p points to rpos[] array; idx of last nonwhitespace char on line
 *   *lalloc_p is >0, and is what the previous read_block call reported
 *   *nlines_p points to the number of lines we saw in the *first* block.
 * On success, returns <eslOK> as above.
 * If the number of lines seen in the block doesn't match the expected number, return <eslEFORMAT>.
 * If no more data remain in the file, return <eslEOF>.
 * If an allocation fails at any point, throw <eslEMEM>. 
 *   
 * Memory for line[][] and llen[] are entirely managed here - not by caller. They are
 * initially allocated on the first block (*lalloc_p == 0); reallocated
 * as needed; and free'd when we reach <EOF> or an error is detected.
 */
static int
read_block(ESL_MSAFILE *afp, char ***line_p, int **llen_p, int **lpos_p, int **rpos_p, int *lalloc_p, int *nlines_p, int *ret_starti)
{
  void  *tmpp;
  char **line;
  int   *llen;
  int   *lpos;
  int   *rpos;
  char  *s;
  int    lalloc;
  int    blen;
  int    nlines = 0;
  int    i;
  int    starti;
  int    status;
  
  afp->errbuf[0] = '\0';
  if (*lalloc_p == 0) 		/* first block? allocate. */
    {
      ESL_ALLOC_CPP(char *, line, sizeof(char *) * 16);
      ESL_ALLOC_CPP(int, llen, sizeof(int)    * 16);
      ESL_ALLOC_CPP(int, lpos, sizeof(int)    * 16);
      ESL_ALLOC_CPP(int, rpos, sizeof(int)    * 16);
      for (i = 0; i < 16; i++) 
	{ line[i] = NULL; llen[i] = 0; lpos[i] = 0; rpos[i] = 0; }
      lalloc = 16;
    }
  else 				/* second or later block? reuse existing arrays   */
    {
      line   = *line_p;
      llen   = *llen_p;
      lpos   = *lpos_p;
      rpos   = *rpos_p;
      lalloc = *lalloc_p;
    }
  
  /* Advance 'til afp->buf contains first line of the block. */
  do { status = msafile_getline(afp); }
  while (status == eslOK && 
	 (is_blankline(afp->buf) || 
	  (*afp->buf == '#' && (strncmp(afp->buf, "#=", 2) != 0))));
  if      (status == eslEOF && *lalloc_p == 0) ESL_XFAIL(eslEOF, afp->errbuf, "parse failed: no alignment data found");
  else if (status != eslOK)                    goto ERROR; /* includes true (normal) EOF, EMEM */

  starti =  afp->linenumber;

  do {
    if (nlines == lalloc) 
      {
	ESL_RALLOC_CPP(char *, line, tmpp, sizeof(char *) * lalloc * 2);
	ESL_RALLOC_CPP(int, llen, tmpp, sizeof(int)    * lalloc * 2);
	ESL_RALLOC_CPP(int, lpos, tmpp, sizeof(int)    * lalloc * 2);
	ESL_RALLOC_CPP(int, rpos, tmpp, sizeof(int)    * lalloc * 2);
	for (i = lalloc; i < lalloc*2; i++) { line[i] = NULL; llen[i] = 0; lpos[i] = 0; rpos[i] = 0; }
	lalloc*=2;
      }

    blen = strlen(afp->buf);
    if (blen > llen[nlines]) ESL_RALLOC_CPP(char, line[nlines], tmpp, sizeof(char) * (blen+1)); /* +1 for \0 */
    strcpy(line[nlines], afp->buf);
    llen[nlines] = blen;

    /* rpos is most efficiently determined in read_block() rather than elsewhere, 
     *  because we know blen here; saves a strlen() elsewhere.
     */
    for (s = line[nlines]+blen-1; isspace(*s); s--) ;
    rpos[nlines] = s - line[nlines];

    nlines++;

    do { status = msafile_getline(afp); }
    while (status == eslOK && (*afp->buf == '#' && strncmp(afp->buf, "#=", 2) != 0)); /* skip comments */
  } while (status == eslOK && ! is_blankline(afp->buf));

  if (status != eslOK && status != eslEOF) goto ERROR; /* EMEM */
  if (*lalloc_p != 0 && *nlines_p != nlines) 
    ESL_XFAIL(eslEFORMAT, afp->errbuf, "parse failed (line %d): expected %d lines in block, saw %d", afp->linenumber, *nlines_p, nlines);

  *line_p     = line;
  *llen_p     = llen; 
  *lpos_p     = lpos;
  *rpos_p     = rpos;
  *lalloc_p   = lalloc;
  *nlines_p   = nlines;
  *ret_starti = starti;
  return eslOK; /* an EOF is turned into OK the first time we see it: so last block read gets dealt with */

 ERROR:  /* includes final EOF, when we try to read another block and fail.  */
  if (line != NULL) { 
    for (i = 0; i < lalloc; i++) 
      if (line[i] != NULL) free(line[i]); 
    free(line);
  }
  if (llen   != NULL) free(llen);
  if (lpos   != NULL) free(lpos);
  if (rpos   != NULL) free(rpos);
  *line_p     = NULL;
  *llen_p     = NULL;
  *lpos_p     = NULL;
  *rpos_p     = NULL;
  *lalloc_p   = 0;
  *nlines_p   = 0;
  *ret_starti = 0;
  return status;	
}

/* First block: determine and store line types, in ltype[0..nlines-1].
 * From that, we know the number of sequences, nseq.
 * From that, we can allocate a new MSA object for <nseq> sequences.
 * Then parse we store all the sequence names in msa->sqname[].
 * This gives us information we will use to validate subsequent blocks,
 * making sure they contain exactly the same line order.
 *
 * We also set lpos[], rpos[] here to the position of the leftmost,
 * rightmost non-whitespace sequence residue character.
 * In the special case of lines with all whitespace data (which SELEX 
 * format allows!), set both lpos[] = -1; we'll catch this
 * as a special case when we need to.
 * 
 * <msa> and <ltype> are allocated here, and must be free'd by caller.
 */
static int
first_selex_block(char *errbuf, int starti, char **line, int *lpos, int *rpos, int nlines, ESL_MSA **ret_msa, int **ret_ltype)
{
  ESL_MSA *msa    = NULL;
  int     *ltype  = NULL;
  int      nseq   = 0;
  int      nrf, ncs, nss, nsa;
  int      has_ss, has_sa;
  int      li, i;
  char    *s, *tok;
  int      n;
  int      status;

  if (errbuf != NULL) errbuf[0] = '\0';

  /* Determine ltype[]; count sequences */
  ESL_ALLOC_CPP(int, ltype, sizeof(int) * nlines);
  nrf = ncs = nss = nsa = 0;
  has_ss = has_sa = FALSE;
  for (nseq = 0, li = 0; li < nlines; li++)
    {
      if      (strncmp(line[li], "#=RF", 4) == 0) { ltype[li] = eslMSA_LINE_RF; nrf++; }
      else if (strncmp(line[li], "#=CS", 4) == 0) { ltype[li] = eslMSA_LINE_CS; ncs++; }
      else if (strncmp(line[li], "#=SS", 4) == 0) { ltype[li] = eslMSA_LINE_SS; nss++;  has_ss = TRUE; }
      else if (strncmp(line[li], "#=SA", 4) == 0) { ltype[li] = eslMSA_LINE_SA; nsa++;  has_sa = TRUE; }
      else                                        { ltype[li] = eslMSA_LINE_SQ; nseq++; nss = nsa = 0; }
      if (nss > 0 && nseq==0) ESL_XFAIL(eslEFORMAT, errbuf, "parse failed (line %d): #=SS must follow a sequence", li+starti);
      if (nsa > 0 && nseq==0) ESL_XFAIL(eslEFORMAT, errbuf, "parse failed (line %d): #=SA must follow a sequence", li+starti);
      if (nrf > 1)            ESL_XFAIL(eslEFORMAT, errbuf, "parse failed (line %d): too many #=RF lines for block", li+starti);
      if (ncs > 1)            ESL_XFAIL(eslEFORMAT, errbuf, "parse failed (line %d): too many #=CS lines for block", li+starti);
      if (nss > 1)            ESL_XFAIL(eslEFORMAT, errbuf, "parse failed (line %d): too many #=SS lines for seq",   li+starti);
      if (nsa > 1)            ESL_XFAIL(eslEFORMAT, errbuf, "parse failed (line %d): too many #=SA lines for seq",   li+starti);
    }

  /* Allocate the MSA, now that we know nseq */
  if ((msa = esl_msa_Create(nseq, -1)) == NULL) { status = eslEMEM; goto ERROR; } 
  if (has_ss) 
    {
      ESL_ALLOC_CPP(char *, msa->ss, sizeof(char *) * nseq); 
      for (i = 0; i < nseq; i++) msa->ss[i] = NULL;
    }
  if (has_sa) 
    {
      ESL_ALLOC_CPP(char *, msa->sa, sizeof(char *) * nseq);
      for (i = 0; i < nseq; i++) msa->sa[i] = NULL;
    }
  msa->nseq = nseq;
  msa->alen = 0;
  /* msa->aseq[], msa->sqname[], msa->ss[], msa->sa[] arrays are all ready (all [i] are NULL) */
  
  for (i = 0, li = 0; li < nlines; li++)
    if (ltype[li] == eslMSA_LINE_SQ)
      {
	s = line[li];
	if (esl_strtok_adv(&s, " \t\n\r", &tok, &n, NULL)     != eslOK) ESL_XEXCEPTION(eslEINCONCEIVABLE, "can't happen");
	if ((status = esl_strdup(tok, -1, &(msa->sqname[i]))) != eslOK) goto ERROR;

	while (*s && isspace(*s)) s++;   /* advance s to first residue */
	lpos[li] = ((*s == '\0') ? -1 : s-line[li]);
	i++;
      }
    else
      {
	for (s = line[li]; *s && !isspace(*s); s++) ; /* advance s past #=XX tag       */
	for (           ;  *s &&  isspace(*s); s++) ; /* advance s to first residue    */
	lpos[li] = ((*s == '\0') ? -1 : s-line[li]);
      }
  *ret_msa   = msa;
  *ret_ltype = ltype;
  return eslOK;

 ERROR:
  if (msa   != NULL) esl_msa_Destroy(msa);
  if (ltype != NULL) free(ltype);
  *ret_msa   = NULL;
  *ret_ltype = NULL;
  return status;
}


/* Subsequent blocks: 
 * validate that lines are coming in same order as first block (including sqname);
 * set lpos[] as in first_selex_block().
 */
static int
other_selex_block(char *errbuf, int starti, char **line, int *lpos, int *rpos, int nlines, ESL_MSA *msa, int *ltype)
{
  char *s, *tok;
  int   i, li;

  /* Compare order of line types. */
  for (li = 0; li < nlines; li++)
    {
      if      (strncmp(line[li], "#=RF", 4) == 0) { if (ltype[li] != eslMSA_LINE_RF) ESL_FAIL(eslEFORMAT, errbuf, "parse failed (line %d): #=RF line isn't in expected order", li+starti); }
      else if (strncmp(line[li], "#=CS", 4) == 0) { if (ltype[li] != eslMSA_LINE_CS) ESL_FAIL(eslEFORMAT, errbuf, "parse failed (line %d): #=CS line isn't in experted order", li+starti); }
      else if (strncmp(line[li], "#=SS", 4) == 0) { if (ltype[li] != eslMSA_LINE_SS) ESL_FAIL(eslEFORMAT, errbuf, "parse failed (line %d): #=SS line isn't in expected order", li+starti); }
      else if (strncmp(line[li], "#=SA", 4) == 0) { if (ltype[li] != eslMSA_LINE_SA) ESL_FAIL(eslEFORMAT, errbuf, "parse failed (line %d): #=SA line isn't in expected order", li+starti); }
      else                                        { if (ltype[li] != eslMSA_LINE_SQ) ESL_FAIL(eslEFORMAT, errbuf, "parse failed (line %d): seq line isn't in expected order",  li+starti); }
    }

  /* Compare order of sequence names, and set lpos[]. */
  for (i = 0, li = 0; li < nlines; li++)
    {
      if (ltype[li] == eslMSA_LINE_SQ)
	{
	  s = line[li];
	  if (esl_strtok(&s, " \t\n\r", &tok) != eslOK) ESL_EXCEPTION(eslEINCONCEIVABLE, "can't happen");
	  if (strcmp(tok, msa->sqname[i])     != 0)     ESL_FAIL(eslEFORMAT, errbuf, "parse failed (line %d): expected seq %s, saw %s",  li+starti, msa->sqname[i], tok);
	  
	  while (*s && isspace(*s)) s++;              /* advance s to first residue */
	  lpos[li] = ((*s == '\0') ? -1 : s-line[li]);
	  i++;
	}
      else
	{
	  for (s = line[li]; *s && !isspace(*s); s++) ; /* advance s past #=XX tag    */
	  for (            ; *s &&  isspace(*s); s++) ; /* advance s to first residue */
	  lpos[li] = ((*s == '\0') ? -1 : s-line[li]);
	}
    }
  return eslOK;
}


static int 
append_selex_block(ESL_MSA *msa, char **line, int *ltype, int *lpos, int *rpos, int nlines)
{
  void *tmpp;
  char *s;
  int   leftmost, rightmost;
  int   li, i;
  int   nleft, ntext, nright;
  int   nadd;
  int   pos;
  int   status;

  /* Determine rightmost, leftmost columns for data */
  /* Watch out for special case of empty data lines: lpos= -1 flag */
  /* Watch out for extra special case where *no* line on block has data! */
  leftmost  = INT_MAX;
  rightmost = -1;
  for (li = 0; li < nlines; li++) {
    leftmost  = (lpos[li] == -1) ? leftmost  : ESL_MIN(leftmost,  lpos[li]);
    rightmost = (lpos[li] == -1) ? rightmost : ESL_MAX(rightmost, rpos[li]);
  }
  if (rightmost == -1) return eslOK; /* extra special case: no data in block at all! */
  nadd = rightmost - leftmost + 1; /* block width in aligned columns */

  for (i = 0, li = 0; li < nlines; li++)
    {
      nleft  = ((lpos[li] != -1) ? lpos[li] - leftmost     : nadd); /* watch special case of all whitespace on data line, lpos>rpos */
      ntext  = ((lpos[li] != -1) ? rpos[li] - lpos[li] + 1 : 0);
      nright = ((lpos[li] != -1) ? rightmost - rpos[li]    : 0);

      if      (ltype[li] == eslMSA_LINE_SQ) { ESL_RALLOC_CPP(char, msa->aseq[i], tmpp, sizeof(char) * (msa->alen + nadd + 1)); s = msa->aseq[i]; i++; }
      else if (ltype[li] == eslMSA_LINE_RF) { ESL_RALLOC_CPP(char, msa->rf,      tmpp, sizeof(char) * (msa->alen + nadd + 1)); s = msa->rf;      }
      else if (ltype[li] == eslMSA_LINE_CS) { ESL_RALLOC_CPP(char, msa->ss_cons, tmpp, sizeof(char) * (msa->alen + nadd + 1)); s = msa->ss_cons; }
      else if (ltype[li] == eslMSA_LINE_SS) { ESL_RALLOC_CPP(char, msa->ss[i-1], tmpp, sizeof(char) * (msa->alen + nadd + 1)); s = msa->ss[i-1]; }
      else if (ltype[li] == eslMSA_LINE_SA) { ESL_RALLOC_CPP(char, msa->sa[i-1], tmpp, sizeof(char) * (msa->alen + nadd + 1)); s = msa->sa[i-1]; }

      for (pos = msa->alen; pos < msa->alen+nleft; pos++) s[pos] = ' ';
      if (ntext > 0) memcpy(s+msa->alen+nleft, line[li]+lpos[li], sizeof(char) * ntext);
      for (pos = msa->alen+nleft+ntext; pos < msa->alen+nadd; pos++) s[pos] = ' ';
      s[pos] = '\0';
    }
  msa->alen += nadd;
  return eslOK;

 ERROR:
  return status;
}
/*------------------- end, selex format -------------------------*/



/*****************************************************************
 * 11. AFA (aligned FASTA) format
 *****************************************************************/

/* read_afa()
 * EPN, Mon Nov  2 17:24:40 2009
 *
 * Purpose:   Parse the one-and-only alignment from an open AFA (aligned
 *            fasta) format alignment file <afp>, leaving the
 *            alignment in <ret_msa>.
 *
 *            The current implementation reads the file one line at a
 *            time. Blank lines are skipped. Lines with '>' as the
 *            first non-whitespace character begin a new sequence,
 *            first word is sequence name, remainder of line is
 *            sequence description. All other lines are sequence lines
 *            currently processed one whitespace-delimited token at a
 *            time (to permit whitespace in the file).  A possibly
 *            more efficient route would be to read each complete
 *            sequence at a time (since AFA is not interleaved) and
 *            write it to the msa in a single step.
 *            
 *            Starting with the second sequence, all sequence lengths
 *            are confirmed to be identical to the length of the first
 *            sequence. If any are not, afp->errbuf is filled, <ret_msa>
 *            is set to NULL and <eslEINVAL> Is returned.
 *
 * Returns:   <eslOK> on success, and the alignment is in <ret_msa>.
 *            If no sequences exist, return <eslEOF> and <ret_msa>
 *            is set to NULL.
 * 
 *            <eslEFORMAT> if parse fails because of a file format
 *            problem, in which case afp->errbuf is set to contain a
 *            formatted message that indicates the cause of the
 *            problem, and <ret_msa> is set to NULL.
 *            
 *            Returns <eslEINVAL> if we're trying to read a digital
 *            alignment, and an invalid residue is found that can't be
 *            digitized.
 * 
 * Throws:    <eslEMEM> on allocation error.
 *
 */
static int
read_afa(ESL_MSAFILE *afp, ESL_MSA **ret_msa)
{
  ESL_MSA   *msa = NULL;
  char      *s;
  int        status;
  int        status2;
  int        seqidx;
  char      *seqname;
  char      *desc;
  char      *text;
  int        len, i;
  char       errbuf2[eslERRBUFSIZE];

  if (feof(afp->f))  { status = eslEOF; goto ERROR; }
  afp->errbuf[0] = '\0';

  /* Initialize allocation of the MSA:
   * make it growable, by giving it an initial blocksize of
   * 16 seqs of 0 length.
   */
#ifdef eslAUGMENT_ALPHABET
  if (afp->do_digital == TRUE && (msa = esl_msa_CreateDigital(afp->abc, 16, -1))  == NULL) 
    { status = eslEMEM; goto ERROR; }

#endif
  if (afp->do_digital == FALSE && (msa = esl_msa_Create(16, -1))  == NULL)
    { status = eslEMEM; goto ERROR; }
  if (msa == NULL)    
    { status = eslEMEM; goto ERROR; }


#ifdef eslAUGMENT_SSI
  /* EPN: not sure if this is appropriate/necessary, we assume only one alignment in AFA files */
  msa->offset = ftello(afp->f);
#endif

  /* Read the alignment file one line at a time.
   */
  while ((status2 = msafile_getline(afp)) == eslOK) 
    {
      s = afp->buf;
      while (*s == ' ' || *s == '\t') s++; /* skip leading whitespace */

      if (*s == '\n' || *s == '\r') continue; /* skip blank lines */
      if (*s == '>') { /* header line */
	/* if nec, make room for the new seq */
	if (msa->nseq >= msa->sqalloc && (status = esl_msa_Expand(msa)) != eslOK) return status; 

	/* store the name (space delimited) */
	s++; /* move past the '>' */
	seqidx = msa->nseq;
	msa->nseq++;
	if (esl_strtok(&s, " \t\n\r", &seqname) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "AFA MSA parse error, problem reading name of sequence %d at line %d\n", seqidx+1, afp->linenumber);
	status = esl_strdup(seqname, -1, &(msa->sqname[seqidx]));

	status = esl_strtok(&s, "\n\r", &desc);
	if     (status == eslOK) status = esl_msa_SetSeqDescription(msa, seqidx, desc);
	else if(status != eslEOL) ESL_XFAIL(eslEFORMAT, afp->errbuf, "AFA MSA parse error, problem reading description of sequence %d at line %d\n", seqidx, afp->linenumber);
	/* else, no description */

	if((seqidx > 1) && (msa->sqlen[(seqidx-1)] != msa->sqlen[0])) { /* make sure the aligned seq we just read is the same length as the first seq we read */
	  ESL_XFAIL(eslEFORMAT, afp->errbuf, "sequence %d length (%ld) is not equal to the expected length (%ld) (the length of first seq in file)", seqidx, msa->sqlen[(seqidx-1)], msa->sqlen[0]);
	}
      }
      else {  /* not a '>' */
	if(msa->nseq == 0) { /* shouldn't happen, we haven't yet seen a '>' */
	  ESL_XFAIL(eslEFORMAT, afp->errbuf, "AFA MSA parse error, first non-whitespace character is not a '>' at line %d\n", afp->linenumber);	  
	}
	/* A sequence line: it doesn't begin with, but may contain, whitespace (' ' or '\t').
	 * We add whitespace-delimited tokens one at a time to the aseq (or ax).
	 * (Note: if we're digitized I think we could use a single call to esl_abc_dsqcat()
	 *  instead of splitting into tokens, which may be slightly more efficient).
	 */
	while(esl_strtok_adv(&s, " \t\n", &text, &len, NULL) == eslOK) 
	  { 
#ifdef eslAUGMENT_ALPHABET
	    if (msa->flags & eslMSA_DIGITAL)
	      {
		if((status = esl_abc_dsqcat(msa->abc, &(msa->ax[seqidx]), &(msa->sqlen[seqidx]), text, len)) != eslOK) {
		  /* invalid char(s), get informative error message */
		  if (esl_abc_ValidateSeq(msa->abc, text, len, afp->errbuf) != eslOK) 
		    ESL_XFAIL(eslEFORMAT, errbuf2, "%s (line %d): %s", msa->sqname[i], afp->linenumber, afp->errbuf);
		}
	      }
#endif
	  if (! (msa->flags & eslMSA_DIGITAL))
	    {
	      status = esl_strcat(&(msa->aseq[seqidx]), msa->sqlen[seqidx], text, len);
	      msa->sqlen[seqidx] += len;
	    } 
	  }
      }
    }

  /* check the length of the final sequence */
  if((msa->nseq > 1) && (msa->sqlen[seqidx] != msa->sqlen[0])) { /* make sure the aligned seq we just read is the same length as the first seq we read */
    ESL_XFAIL(eslEINVAL, afp->errbuf, "sequence %d length (%ld) is not equal to the expected length (%ld) (the length of first seq in file)", seqidx+1, msa->sqlen[seqidx], msa->sqlen[0]);
  }

  if (status2 == eslEMEM) ESL_XFAIL(status, afp->errbuf, "out of memory");
  if (status2 != eslEOF)  ESL_XFAIL(status, afp->errbuf, "unexpected error reading AFA alignment");

  /* Verify the msa */
  if (verify_parse(msa, afp->errbuf) != eslOK) { status = eslEFORMAT; goto ERROR; } 

  /* if alignment is empty set <ret_msa> to NULL and return eslEOF, (verification still works in this case) */
  if (msa->nseq == 0) { status = eslEOF; goto ERROR; }

  if (ret_msa != NULL) *ret_msa = msa; else esl_msa_Destroy(msa);
  return eslOK;

 ERROR:
  if (msa != NULL)      esl_msa_Destroy(msa);
  if (ret_msa != NULL) *ret_msa = NULL;
  return status;
}

/*---------------------- end, AFA format ------------------------*/

/*****************************************************************
 * 12.5. galosh profile format (from profilic)
 *****************************************************************/

/* profillic_read_profile()
 * Paul T Edlefsen   paul@galosh.org   June 29, 2011.
 *
 * Purpose: Parse the Profile HMM from an open galosh profile format
 *            file (from profillic) <afp>, leaving the profile in
 *            <ret_profile>.
 *
 * Returns:   <eslOK> on success, and the alignment is in <ret_profile>.
 * 
 *            <eslEFORMAT> if parse fails because of a file format
 *            problem, in which case afp->errbuf is set to contain a
 *            formatted message that indicates the cause of the
 *            problem, and <ret_profile> is unaffected.
 *            
 * Throws:    <eslEMEM> on allocation error.
 *
 */
template <typename ProfileType>
static int
profillic_read_profile(ESL_MSAFILE *afp, ESL_MSA **ret_msa, ProfileType * profile_ptr)
{
  ESL_MSA   *msa = NULL;
  char *buf;
  long len;
  string profile_string;
  int        seqidx;
  int status;
  char       errbuf2[eslERRBUFSIZE];

  const char * const seqname = "Galosh Profile Consensus";
  const char * const msaname = "Galosh Profile";
  uint32_t profile_length;
  galosh::Sequence<typename ProfileType::ProfileResidueType> consensus_sequence;
  stringstream tmp_consensus_output_stream;

  uint32_t pos_i;

  if (profile_ptr == NULL)  { ESL_EXCEPTION(eslEINCONCEIVABLE, "profile_ptr is NULL in profillic_read_profile(..)!"); }
  if (feof(afp->f))  { status = eslEOF; goto ERROR; }
  afp->errbuf[0] = '\0';

  // Read in the galosh profile (from profillic)
  fseek( afp->f, 0, SEEK_END ); // go to the end
  len = ftell( afp->f ); // get the position at the end (length)
  fseek( afp->f, 0, SEEK_SET ); // go to the beginning again.

  ESL_ALLOC_CPP( char, buf, sizeof( char ) * len ); //malloc buffer
  fread( buf, len, 1, afp->f ); //read into buffer

  profile_string = buf;
  profile_ptr->fromString( profile_string );
  if (buf)      free(buf);
  // TODO: WHY WON'T THIS WORK?  See HACKs in profillic-hmmbuild.cpp to work around it.
  fseek( afp->f, 0, SEEK_END ); // go to the end (to signal there's no more profiles in the file, the next time we come to this function)

  // Calculate the consensus sequence.
  profile_length = profile_ptr->length();
  consensus_sequence.reinitialize( profile_length );
  for( pos_i = 0; pos_i < profile_length; pos_i++ ) {
    consensus_sequence[ pos_i ] =
      ( *profile_ptr )[ pos_i ][ galosh::Emission::Match ].maximumValueType();
  }
  tmp_consensus_output_stream << consensus_sequence;

  /* Initialize allocation of the MSA:
   * make it growable, by giving it an initial blocksize of
   * 16 seqs of 0 length.
   */
#ifdef eslAUGMENT_ALPHABET
  if (afp->do_digital == TRUE && (msa = esl_msa_CreateDigital(afp->abc, 16, -1))  == NULL) 
    { status = eslEMEM; goto ERROR; }

#endif
  if (afp->do_digital == FALSE && (msa = esl_msa_Create(16, -1))  == NULL)
    { status = eslEMEM; goto ERROR; }
  if (msa == NULL)    
    { status = eslEMEM; goto ERROR; }


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
      if((status = esl_abc_dsqcat(msa->abc, &(msa->ax[seqidx]), &(msa->sqlen[seqidx]), tmp_consensus_output_stream.str().c_str(), profile_length)) != eslOK) {
        /* invalid char(s), get informative error message */
        if (esl_abc_ValidateSeq(msa->abc, tmp_consensus_output_stream.str().c_str(), profile_length, afp->errbuf) != eslOK) 
          ESL_XFAIL(eslEFORMAT, errbuf2, "%s (line %d): %s", msa->sqname[0], afp->linenumber, afp->errbuf);
      }
    }
#endif
  if (! (msa->flags & eslMSA_DIGITAL))
    {
      status = esl_strcat(&(msa->aseq[seqidx]), 0, tmp_consensus_output_stream.str().c_str(), profile_length);
      msa->sqlen[seqidx] = profile_length;
    } 

  /// .... OR read in a fasta file of sequences too.
  // TODO: (Optional?) Set msa->name to the name of the profile (file?)
  esl_strdup(msaname, -1, &(msa->name));
  // TODO: make sure eslMSA_HASWGTS is FALSE .. OR set it to TRUE and set msa->wgt[idx] to 1.0.
  // NOTE: Could have secondary structure (per sequence) too. msa->ss[0]. msa->sslen[0] should be the same as msa->sqlen[0].
  // TODO: Investigate what msa->sa and msa->pp are for.

  /* Give the newly parsed MSA a good
   * going-over, and finalize the fields of the MSA data structure.
   * verify_parse will fill in errbuf if it sees a problem.
   */
  if (verify_parse(msa, afp->errbuf) != eslOK) { status = eslEFORMAT; goto ERROR; } 

  if (ret_msa != NULL) *ret_msa = msa; else esl_msa_Destroy(msa);
  return eslOK;

 ERROR:
  if (msa != NULL)      esl_msa_Destroy(msa);
  if (ret_msa != NULL) *ret_msa = NULL;
  return status;
}

/*---------------------- end, galosh profile format (from profillic)-------*/


/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.0; March 2010
 * Copyright (C) 2010 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/

#endif // __GALOSH_PROFILLICESLMSA_HPP__
