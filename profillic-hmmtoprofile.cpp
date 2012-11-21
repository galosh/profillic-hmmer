/**
 * \file profillic-hmmtoprofile.cpp
 * \brief
 * Convert HMM to galosh profile
 * \details
<pre>
# profillic-hmmtoprofile :: convert HMM to galosh profile
# profillic-hmmer 1.0a (July 2011); http://galosh.org/
# Copyright (C) 2011 Paul T. Edlefsen, Fred Hutchinson Cancer Research Center.
# HMMER 3.1dev (November 2011); http://hmmer.org/
# Copyright (C) 2011 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Usage: profillic-hmmtoprofile [-options] <input hmmfile> <output galosh profile>

Options:
  -h : show brief help on version and usage
</pre>
 */
extern "C" {
#include "p7_config.h"
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

extern "C" {
#include "easel.h"
#include "esl_getopts.h"
/// \note TAH 8/12 workaround to avoid C++ keyword "new" in esl_msa.h
#define new _new
#include "hmmer.h"
#undef new
}

/* ////////////// For profillic-hmmer ////////////////////////////////// */
#include "profillic-hmmer.hpp"

#include <iostream>

// Updated notices:
#define PROFILLIC_HMMER_VERSION "1.0a"
#define PROFILLIC_HMMER_DATE "July 2011"
#define PROFILLIC_HMMER_COPYRIGHT "Copyright (C) 2011 Paul T. Edlefsen, Fred Hutchinson Cancer Research Center."
#define PROFILLIC_HMMER_URL "http://galosh.org/"

// Modified from hmmer.c p7_banner(..):
/* Version info - set once for whole package in configure.ac
 */
/*****************************************************************
 * 1. Miscellaneous functions for H3
 *****************************************************************/

/**
 * <pre>
 * Function:  p7_banner()
 * Synopsis:  print standard HMMER application output header
 * Incept:    SRE, Wed May 23 10:45:53 2007 [Janelia]
 *
 * Purpose:   Print the standard HMMER command line application banner
 *            to <fp>, constructing it from <progname> (the name of the
 *            program) and a short one-line description <banner>.
 *            For example, 
 *            <p7_banner(stdout, "hmmsim", "collect profile HMM score distributions");>
 *            might result in:
 *            
 *            \begin{cchunk}
 *            # hmmsim :: collect profile HMM score distributions
 *            # HMMER 3.0 (May 2007)
 *            # Copyright (C) 2004-2007 HHMI Janelia Farm Research Campus
 *            # Freely licensed under the Janelia Software License.
 *            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *            \end{cchunk}
 *              
 *            <progname> would typically be an application's
 *            <argv[0]>, rather than a fixed string. This allows the
 *            program to be renamed, or called under different names
 *            via symlinks. Any path in the <progname> is discarded;
 *            for instance, if <progname> is "/usr/local/bin/hmmsim",
 *            "hmmsim" is used as the program name.
 *            
 * Note:    
 *    Needs to pick up preprocessor #define's from p7_config.h,
 *    as set by ./configure:
 *            
 *    symbol          example
 *    ------          ----------------
 *    HMMER_VERSION   "3.0"
 *    HMMER_DATE      "May 2007"
 *    HMMER_COPYRIGHT "Copyright (C) 2004-2007 HHMI Janelia Farm Research Campus"
 *    HMMER_LICENSE   "Freely licensed under the Janelia Software License."
 *
 * Returns:   (void)
 * </pre>
 */
void
profillic_p7_banner(FILE *fp, char *progname, char *banner)
{
  char *appname = NULL;

  if (esl_FileTail(progname, FALSE, &appname) != eslOK) appname = progname;

  fprintf(fp, "# %s :: %s\n", appname, banner);
  fprintf(fp, "# profillic-hmmer %s (%s); %s\n", PROFILLIC_HMMER_VERSION, PROFILLIC_HMMER_DATE, PROFILLIC_HMMER_URL);
  fprintf(fp, "# %s\n", PROFILLIC_HMMER_COPYRIGHT);
  fprintf(fp, "# HMMER %s (%s); %s\n", HMMER_VERSION, HMMER_DATE, HMMER_URL);
  fprintf(fp, "# %s\n", HMMER_COPYRIGHT);
  fprintf(fp, "# %s\n", HMMER_LICENSE);
  fprintf(fp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

  if (appname != NULL) free(appname);
  return;
}

template <typename ProfileType>
int
convert_to_galosh_profile ( P7_HMM * hmm, ProfileType & profile )
{
  typedef typename galosh::profile_traits<ProfileType>::ResidueType ResidueType;

  int status;

  uint32_t pos_i; // Position in profile.  Corresponds to one less than match state pos in HMM.
  uint32_t res_i;
  ESL_DSQ hmmer_digitized_residue;

  /* How many match states in the HMM? */
  if( hmm->M == 0 ) { status = eslENORESULT; goto ERROR; }
  profile.reinitialize( static_cast<uint32_t>( hmm->M ) );

  profile.zero();

  /// \note NOTE that HMMER3 has a slightly different model, starting in
  /// Begin rather than in preAlign, and with 3 legal transitions out
  /// of Begin (one of these is to PreAlign).  The galosh profile model
  /// begins in preAlign and transitions to Begin, and from there to
  /// either Match or Delete.  One implication is that galosh profiles
  /// enforce t[ 0 ][ p7H_MI ] to be the same as t[ 0 ][ p7H_II ], but
  /// HMMER3 does not.  Another way to say this is that H3 uses affine
  /// pre-aligns, and prohibits pre-align -to- delete transitions,
  /// whereas galosh / profillic uses non-affine pre-aligns and allows
  /// pre-align->delete.

  // fromPreAlign
  profile[ galosh::Transition::fromPreAlign ][ galosh::TransitionFromPreAlign::toPreAlign ] =
    hmm->t[ 0 ][ p7H_II ];
  profile[ galosh::Transition::fromPreAlign ][ galosh::TransitionFromPreAlign::toBegin ] =
    hmm->t[ 0 ][ p7H_IM ];
  for( res_i = 0; res_i < seqan::ValueSize<ResidueType>::VALUE; res_i++ ) {
    hmmer_digitized_residue =
      esl_abc_DigitizeSymbol( hmm->abc, static_cast<char>( ResidueType( res_i ) ) );
    // See below where it says "TODO/NOTE"..
    profile[ galosh::Emission::PreAlignInsertion ][ res_i ] =
      hmm->ins[ 0 ][ hmmer_digitized_residue ];
  }

  // fromBegin
  profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toMatch ] =
    ( hmm->t[ 0 ][ p7H_MM ] / ( 1.0 - hmm->t[ 0 ][ p7H_MI ] ) );
  profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toDeletion ] =
    ( 1.0 - profile[ galosh::Transition::fromBegin ][ galosh::TransitionFromBegin::toMatch ] );

  for( pos_i = 0; pos_i < profile.length(); pos_i++ ) {
//    if( be_verbose ) {
//      cout << '.';
//      cout.flush();
//    }
    // TODO: If this is too slow, memoize the ResidueType( res_i )s.
    for( res_i = 0; res_i < seqan::ValueSize<ResidueType>::VALUE; res_i++ ) {
      hmmer_digitized_residue =
        esl_abc_DigitizeSymbol( hmm->abc, static_cast<char>( ResidueType( res_i ) ) );
      profile[ pos_i ][ galosh::Emission::Match ][ res_i ] =
        hmm->mat[ pos_i + 1 ][ hmmer_digitized_residue ];
      if( pos_i == ( profile.length() - 1 ) ) {
        // Use post-align insertions
        profile[ galosh::Emission::PostAlignInsertion ][ res_i ] =
          hmm->ins[ pos_i + 1 ][ hmmer_digitized_residue ];
      } else { // if this is the last position (use post-align insertions) .. else ..
        profile[ galosh::Emission::Insertion ][ res_i ] +=
          hmm->ins[ pos_i + 1 ][ hmmer_digitized_residue ];
      } // End if this is the last position (use post-align insertions) .. else ..
    } // End foreach res_i
    if( pos_i == ( profile.length() - 1 ) ) {
      // Use post-align insertions
      profile[ galosh::Transition::fromPostAlign ][ galosh::TransitionFromPostAlign::toTerminal ] =
        hmm->t[ pos_i + 1 ][ p7H_IM ];
      profile[ galosh::Transition::fromPostAlign ][ galosh::TransitionFromPostAlign::toPostAlign ] =
        ( 1.0 - profile[ galosh::Transition::fromPostAlign ][ galosh::TransitionFromPostAlign::toTerminal ] );
    } else {  // if this is the last position (use post-align insertions) .. else ..
      profile[ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toMatch ] +=
        hmm->t[ pos_i + 1 ][ p7H_MM ];
      profile[ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toInsertion ] +=
        hmm->t[ pos_i + 1 ][ p7H_MI ];
      profile[ galosh::Transition::fromMatch ][ galosh::TransitionFromMatch::toDeletion ] +=
        hmm->t[ pos_i + 1 ][ p7H_MD ];
  
      profile[ galosh::Transition::fromInsertion ][ galosh::TransitionFromInsertion::toMatch ] +=
        hmm->t[ pos_i + 1 ][ p7H_IM ];
      profile[ galosh::Transition::fromInsertion ][ galosh::TransitionFromInsertion::toInsertion ] +=
        hmm->t[ pos_i + 1 ][ p7H_II ];
      profile[ galosh::Transition::fromDeletion ][ galosh::TransitionFromDeletion::toMatch ] +=
        hmm->t[ pos_i + 1 ][ p7H_DM ];
      profile[ galosh::Transition::fromDeletion ][ galosh::TransitionFromDeletion::toDeletion ] +=
        hmm->t[ pos_i + 1 ][ p7H_DD ];
    } // End if this is the last position (use post-align insertions) .. else ..
  } // End foreach pos_i

  // Normalize with 0 as the minimum value we'll allow.  Note that in
  // profillic and profuse, it's generally 1E-5, so when the profile
  // is read in by those programs, it might be slightly altered.
  profile.normalize( 0 );
  return eslOK;

 ERROR:
  return status;
} // convert_to_galosh_profile (..)

/* ////////////// End profillic-hmmer ////////////////////////////////// */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles   reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL, NULL, "show brief help on version and usage",            0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <input hmmfile> <output galosh profile>";
static char banner[] = "convert HMM to galosh profile";
/**
 * \fn int main(int argc,char **argv)
 * main driver
 *
 */
int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go	   = NULL;      /**< command line processing                   */
  ESL_ALPHABET    *abc     = NULL;
  char            *hmmfile = NULL;
  char            *outhmmfile = NULL;
  P7_HMMFILE      *hfp     = NULL;
  P7_HMM          *hmm     = NULL;
  P7_BG           *bg      = NULL;
  int              nhmm;	
  double           x;
  float            KL;
  int              status;
  char             errbuf[eslERRBUFSIZE];

  char        errmsg[eslERRBUFSIZE];

  /* Process the command line options.
   */
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK || 
      esl_opt_VerifyConfig(go)               != eslOK)
    {
      printf("Failed to parse command line: %s\n", go->errbuf);
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      profillic_p7_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nOptions:");
      esl_opt_DisplayHelp(stdout, go, 0, 2, 80); /* 0=docgroup, 2 = indentation; 80=textwidth*/
      exit(0);
    }
  if (esl_opt_ArgNumber(go) != 2) 
    {
      puts("Incorrect number of command line arguments.");
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  if ((hmmfile = esl_opt_GetArg(go, 1)) == NULL) 
    {
      puts("Failed to read <input hmmfile> argument from command line.");
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  if ((outhmmfile = esl_opt_GetArg(go, 2)) == NULL) 
    {
      puts("Failed to read <output hmmfile> argument from command line.");
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  profillic_p7_banner(stdout, argv[0], banner);
  
  /* Initializations: open the input HMM file for reading
   */
  status = p7_hmmfile_OpenE(hmmfile, NULL, &hfp, errbuf);
  if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", hmmfile, errbuf);
  else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                hmmfile, errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, hmmfile, errbuf);  

  /** Main body: read HMMs one at a time, print one line of stats
   */
  printf("#\n");
  printf("# %-4s %-20s %-12s %8s %8s %6s %6s %6s %6s %6s\n", "idx",  "name",                 "accession",    "nseq",     "eff_nseq", "M",      "relent", "info",   "p relE", "compKL");
  printf("# %-4s %-20s %-12s %8s %8s %6s %6s %6s %6s %6s\n", "----", "--------------------", "------------", "--------", "--------", "------", "------", "------", "------", "------");

  nhmm = 0;
  if ((status = p7_hmmfile_Read(hfp, &abc, &hmm)) != eslEOF) 
    {
      if      (status == eslEOD)       esl_fatal("read failed, HMM file %s may be truncated?", hmmfile);
      else if (status == eslEFORMAT)   esl_fatal("bad file format in HMM file %s",             hmmfile);
      else if (status == eslEINCOMPAT) esl_fatal("HMM file %s contains different alphabets",   hmmfile);
      else if (status != eslOK)        esl_fatal("Unexpected error in reading HMMs from %s",   hmmfile);
      nhmm++;

      if (bg == NULL) bg = p7_bg_Create(abc);

      if( abc->type == eslDNA ) {
        galosh::ProfileTreeRoot<seqan::Dna, floatrealspace> profile;
        if( (status = convert_to_galosh_profile( hmm, profile )) != eslOK ) esl_fatal("Unexpected error in converting HMM from file %s to a dna galosh profile",   hmmfile);
        std::ofstream fs ( outhmmfile );
    
        if( !fs.is_open() ) {
          esl_fatal("Unexpected error in opening the file %s for writing", outhmmfile);
        } else {
          fs << profile;
          fs.close();
        }
      } else if( abc->type == eslAMINO ) {
        galosh::ProfileTreeRoot<seqan::AminoAcid20, floatrealspace> profile;
        if( (status = convert_to_galosh_profile( hmm, profile )) != eslOK ) esl_fatal("Unexpected error in converting HMM from file %s to an amino galosh profile",   hmmfile);
        std::ofstream fs ( outhmmfile );
    
        if( !fs.is_open() ) {
          esl_fatal("Unexpected error in opening the file %s for writing", outhmmfile);
        } else {
          fs << profile;
          fs.close();
        }
      } else {
        ESL_EXCEPTION(eslEUNIMPLEMENTED, "Sorry, at present the profillic-hmmtoprofile software can only handle amino and dna.");
      }
  
      p7_MeanPositionRelativeEntropy(hmm, bg, &x); 
      p7_hmm_CompositionKLDist(hmm, bg, &KL, NULL);

      printf("%-6d %-20s %-12s %8d %8.2f %6d %6.2f %6.2f %6.2f %6.2f\n",
	     nhmm,
	     hmm->name,
	     hmm->acc == NULL ? "-" : hmm->acc,
	     hmm->nseq,
	     hmm->eff_nseq,
	     hmm->M,
	     p7_MeanMatchRelativeEntropy(hmm, bg),
	     p7_MeanMatchInfo(hmm, bg),
	     x,
	     KL);

	     /*	     p7_MeanForwardScore(hmm, bg)); */

      p7_hmm_Destroy(hmm);
    }

  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  p7_hmmfile_Close(hfp);
  esl_getopts_Destroy(go);
  exit(0);
}
