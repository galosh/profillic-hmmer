/**
 * \file profillic-hmmcalibrate.cpp
 * \brief Calibrate HMM search statistics 
 * \details
 * <pre>
# profillic-hmmcalibrate :: calibrate HMM search statistics
# profillic-hmmer 1.0a (July 2011); http://galosh.org/
# Copyright (C) 2011 Paul T. Edlefsen, Fred Hutchinson Cancer Research Center.
# HMMER 3.1dev (November 2011); http://hmmer.org/
# Copyright (C) 2011 Howard Hughes Medical Institute.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Usage: profillic-hmmcalibrate [-options] <input hmmfile> <output hmmfile>

Options:
  -h         : show brief help on version and usage
  --seed <n> : set RNG seed to <n> (if 0: one-time arbitrary seed)  [42]  (n>=0)
 * </pre>
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
  /// \note TAH 8/12 Workaround to avoid use of C++ keyword "new" in esl_msa.h
#define new _new
#include "hmmer.h"
#undef new
}

/* ////////////// For profillic-hmmer ///////////////////////////////// */
#include "profillic-hmmer.hpp"
//#include "profillic-p7_builder.hpp"

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
/* ////////////// End profillic-hmmer ////////////////////////////////// */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles   reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL, NULL, "show brief help on version and usage",            0 },
  { "--seed",     eslARG_INT,   "42", NULL, "n>=0",     NULL,      NULL,    NULL, "set RNG seed to <n> (if 0: one-time arbitrary seed)",   8 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <input hmmfile> <output hmmfile>";
static char banner[] = "calibrate HMM search statistics";
/**
 * int main(int argc, char **argv) 
 * main driver
 *
 */
int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go	   = NULL;      /* command line processing                   */
  ESL_ALPHABET    *abc     = NULL;
  char            *hmmfile = NULL;
  char            *outhmmfile = NULL;
  P7_HMMFILE      *hfp     = NULL;
  FILE         *outhmmfp;          /* HMM output file handle                  */
  P7_HMM          *hmm     = NULL;
  P7_BG           *bg      = NULL;
  int              nhmm;	
  double           x;
  float            KL;
  int              status;
  char             errbuf[eslERRBUFSIZE];

  char        errmsg[eslERRBUFSIZE];

  /* Run-to-run variation due to random number generation                                          */
  int         seed;
  ESL_RANDOMNESS      *r;	         /* RNG for E-value calibration simulations                */
  int                  do_reseeding;	 /* TRUE to reseed, making results reproducible            */

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
  if (esl_opt_IsUsed(go, "--seed"))  {
    if (esl_opt_GetInteger(go, "--seed") == 0) printf("# random number seed:               one-time arbitrary\n");
    else                                       printf("# random number seed set to:        %d\n", esl_opt_GetInteger(go, "--seed"));
  }
  
  /* Initializations: open the input HMM file for reading
   */
  status = p7_hmmfile_OpenE(hmmfile, NULL, &hfp, errbuf);
  if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", hmmfile, errbuf);
  else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                hmmfile, errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, hmmfile, errbuf);  

  /* Initializations: open the output HMM file for writing
   */
  if ((outhmmfp = fopen(outhmmfile, "w")) == NULL) ESL_FAIL(status, errmsg, "Failed to open HMM file %s for writing", outhmmfile);

  /* Normally we reinitialize the RNG to original seed before calibrating each model.
   * This eliminates run-to-run variation.
   * As a special case, seed==0 means choose an arbitrary seed and shut off the
   * reinitialization; this allows run-to-run variation.
   */
  seed       = esl_opt_GetInteger(go, "--seed");
  r            = esl_randomness_CreateFast(seed);
  do_reseeding = (seed == 0) ? FALSE : TRUE;

  /* Main body: read HMMs one at a time, print one line of stats
   */
  printf("#\n");
  printf("# %-4s %-20s %-12s %8s %8s %6s %6s %6s %6s %6s\n", "idx",  "name",                 "accession",    "nseq",     "eff_nseq", "M",      "relent", "info",   "p relE", "compKL");
  printf("# %-4s %-20s %-12s %8s %8s %6s %6s %6s %6s %6s\n", "----", "--------------------", "------------", "--------", "--------", "------", "------", "------", "------", "------");

  nhmm = 0;
  while ((status = p7_hmmfile_Read(hfp, &abc, &hmm)) != eslEOF) 
    {
      if      (status == eslEOD)       esl_fatal("read failed, HMM file %s may be truncated?", hmmfile);
      else if (status == eslEFORMAT)   esl_fatal("bad file format in HMM file %s",             hmmfile);
      else if (status == eslEINCOMPAT) esl_fatal("HMM file %s contains different alphabets",   hmmfile);
      else if (status != eslOK)        esl_fatal("Unexpected error in reading HMMs from %s",   hmmfile);
      nhmm++;

      if (bg == NULL) bg = p7_bg_Create(abc);

      /// \todo Add use of profillic-p7_builder and command-line args to control calibration.
      if ((status = p7_Calibrate(hmm, NULL, &r, &bg, NULL, NULL)) != eslOK) esl_fatal("Unexpected error in calibrating the hmm");

      if( do_reseeding ) {
        // For next time, reset the RNG to what it was this time..
        esl_randomness_Init(r, esl_randomness_GetSeed(r));
      }

      if ((status = p7_hmm_Validate(hmm, errmsg, 0.0001))       != eslOK) return status;
      if ((status = p7_hmmfile_WriteASCII(outhmmfp, -1, hmm)) != eslOK) ESL_FAIL(status, errmsg, "HMM save failed");
  
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
  if (outhmmfp != NULL) fclose(outhmmfp);
 esl_getopts_Destroy(go);
  exit(0);
}
