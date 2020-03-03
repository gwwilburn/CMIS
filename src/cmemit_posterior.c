#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_vectorops.h"

#include "p7_config.h"
#include "hmmer.h"

#include "infernal.h"
#include "config.h"

/* declaration of internal functions */
int generate_msa(CM_t *cm, ESL_SQ **sq, ESL_RANDOMNESS *rng, int nseq, int R, ESL_MSA **ret_msa, char *errbuf);

static ESL_OPTIONS options[] = {
        /* name       type        default    env range togs  reqs  incomp            help                                             docgroup */
        { "-h",       eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL,            "help; show brief info on version and usage",   1 },

        /* Options forcing which alphabet we're working in (normally autodetected) */
        { "--amino",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, "--dna,--rna",   "We are dealing with protein sequences",        2 },
        { "--dna",    eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, "--amino,--rna", "We are dealing with dna sequences",            2 },
        { "--rna",    eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, "--amino,--dna", "We are dealing with rna sequences",            2 },

        /* options for controlling sampling */
        { "-R",       eslARG_INT,     "10",  NULL, NULL, NULL, NULL, NULL,            "Number of CM paths sampled per sequence",      3 },
        { "-s",       eslARG_INT,      "0",   NULL, NULL, NULL, NULL, NULL,           "Set random number seed to <n>",                3 },

        /* control of output */

        { 0,0,0,0,0,0,0,0,0,0 },
};

static char banner[] = "Emit alignments from P_cm (path | sequence)";
static char usage[]  = "[-options] <cmfile> <seqfile> <msa_outfile>";

static void
cmdline_failure(char *argv0, char *format, ...)
{
        va_list argp;
        printf("\nERROR: ");
        va_start(argp, format);
        vfprintf(stderr, format, argp);
        va_end(argp);
        esl_usage(stdout, argv0, usage);
        printf("\nTo see more help on available options, do %s -h\n\n", argv0);
        exit(1);
}

int main (int argc, char *argv[]) {
   ESL_GETOPTS      *go               = NULL;                 /* application configuration                */
   ESL_ALPHABET     *abc              = NULL;                 /* biological alphabet                      */
   ESL_RANDOMNESS   *rng              = NULL;                 /* random number generator                  */
   char             *cmfile           = NULL;                 /* input cm filepath                        */
   char             *seqfile          = NULL;                 /* input seq filepath                       */
   char             *msafile          = NULL;                 /* output msa filepath                      */
   CM_FILE          *cmfp             = NULL;                 /* open input CM file stream                */
   ESL_SQFILE       *sqfp             = NULL;                 /* open seq file stream                     */
   CM_t             *cm               = NULL;                 /* cm                                       */
   ESL_SQ          **sq               = NULL;                 /* array of sequences                       */
   int               format           = eslSQFILE_UNKNOWN;    /* seq file format                          */
   FILE             *afp              = NULL;                 /* output alignment file (-A)               */
   ESL_MSA          *msa              = NULL;                 /* alignment output object                  */
   int               outfmt           = eslMSAFILE_STOCKHOLM; /* alignment output format                  */
   int               n                = 0;                    /* sequence index                           */
   int               nseq             = 0;                    /* total number of seqs in seq file         */
   int               R;                                       /* number of paths sampled per sequence     */
   int               status;                                  /* easel return code                        */
   char              errbuf[eslERRBUFSIZE];


   /* parse command line */
   go = esl_getopts_Create(options);
   if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) cmdline_failure(argv[0], "Failed to parse command line: %s\n", go->errbuf);
   if (esl_opt_VerifyConfig(go)               != eslOK) cmdline_failure(argv[0], "Error in app configuration:   %s\n", go->errbuf);

   if (esl_opt_GetBoolean(go, "-h") ) {
      esl_banner(stdout, argv[0], banner);
      esl_usage (stdout, argv[0], usage);
      puts("\n where options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      puts("\n Alphabet options:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80);
      puts("\n Sampling options:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);

      exit(0);
   }

   /* read arguments */
   if (esl_opt_ArgNumber(go) != 3)  cmdline_failure(argv[0], "Incorrect number of command line arguments.\n", go->errbuf);

   cmfile    =  esl_opt_GetArg(go, 1);
   seqfile   =  esl_opt_GetArg(go, 2);
   msafile   =  esl_opt_GetArg(go, 3);


   /* if user has defined an alphabet we define it here */
   if        (esl_opt_GetBoolean(go, "--amino"))       abc = esl_alphabet_Create(eslAMINO);
   else if   (esl_opt_GetBoolean(go, "--rna"))         abc = esl_alphabet_Create(eslRNA);
   else if   (esl_opt_GetBoolean(go, "--dna"))         abc = esl_alphabet_Create(eslDNA);

   /* get number of paths */
   R = esl_opt_GetInteger(go, "-R");

   /* create random number generator */
   rng = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

   /* open the .cm file */
   status = cm_file_Open(cmfile, NULL, FALSE, &cmfp, errbuf);
   if      (status == eslENOTFOUND) cm_Fail("File existence/permissions problem in trying to open CM file %s.\n%s\n", cmfile, errbuf);
   else if (status == eslEFORMAT)   cm_Fail("File format problem in trying to open CM file %s.\n%s\n",                cmfile, errbuf);
   else if (status != eslOK)        cm_Fail("Unexpected error %d in opening CM file %s.\n%s\n",               status, cmfile, errbuf);

   /* read first cm in file */
   status = cm_file_Read(cmfp, TRUE, &abc, &cm);
   if      (status == eslEFORMAT)   cm_Fail("Bad file format in CM file %s:\n%s\n",          cmfp->fname, cmfp->errbuf);
   else if (status == eslEINCOMPAT) cm_Fail("CM in %s is not in the expected %s alphabet\n", cmfp->fname, esl_abc_DecodeType(abc->type));
   else if (status == eslEOF)       cm_Fail("Empty CM file %s? No CM data found.\n",         cmfp->fname);
   else if (status != eslOK)        cm_Fail("Unexpected error in reading CMs from %s\n",     cmfp->fname);

   cm_file_Close(cmfp);

   /* configure CM */
   cm->config_opts |= CM_CONFIG_SCANMX;
   cm_Configure(cm, errbuf, -1);

   /* set search options */
   cm->search_opts |= CM_SEARCH_INSIDE;
   cm->search_opts |= CM_SEARCH_NONBANDED;

   /* initlalize logsum lookup table magic */
   init_ilogsum();
   FLogsumInit();

   /* print all configuration flags */
   DumpCMFlags(stdout, cm);

   /* open the sequence file */
   status = esl_sqfile_OpenDigital(abc, seqfile, format, NULL, &sqfp);
   if      (status == eslENOTFOUND) p7_Fail("No such file.");
   else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
   else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
   else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);

   /* read sequences into array */
   ESL_REALLOC(sq, sizeof(ESL_SQ *) * (n + 1));
   sq[n] = esl_sq_CreateDigital(abc);
   while ((status = esl_sqio_Read(sqfp, sq[n+nseq])) == eslOK) {
      nseq++;
      ESL_REALLOC(sq, sizeof(ESL_SQ *) * (n+nseq+1));
      sq[n+nseq] = esl_sq_CreateDigital(abc);
   }

   /* error handling and cleanup */
   if      (status == eslEFORMAT) esl_fatal("Parse failed (sequence file %s):\n%s\n",
                                          sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
   else if (status != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s",
                                            status, sqfp->filename);

   if ((status = generate_msa(cm, sq, rng, nseq, R, &msa, errbuf) != eslOK))  p7_Fail("issue running generate_msa(), returned code %d", status);

   /* write MSA to output file */
   if ( (afp = fopen(msafile,"w"))== NULL)  p7_Fail("Failed to open alignment file %s for writing\n", msafile);
   esl_msafile_Write(afp, msa, outfmt);
   fclose(afp);

   /* clean up and return */
   for (n = 0; n < nseq+1; n++)
   {
      esl_sq_Destroy(sq[n]);
   }

   esl_sqfile_Close(sqfp);
   esl_msa_Destroy(msa);
   free(sq);
   FreeCM(cm);
   esl_randomness_Destroy(rng);
   esl_alphabet_Destroy(abc);
   esl_getopts_Destroy(go);
   return 0;

   ERROR:
      return status;
}

int generate_msa(CM_t *cm, ESL_SQ **sq, ESL_RANDOMNESS *rng, int nseq, int R, ESL_MSA **ret_msa, char *errbuf) {
   int           nrow = R * nseq;             /* total number of paths we will generate        */
   int           n;                           /* sequence index [0,...,nseq-1]                 */
   int           r;                           /* samp;e index [0,...,R-1]                      */
   int           idx;                         /* combined sequence-sample index [0,...,nrow-1] */
   float         insc;                        /* inside score                                  */
   float         sc;                          /* S_cm (x, pi)                                  */
   ESL_MSA      *msa = NULL;                  /* msa we're building                            */
   Parsetree_t **pstr = NULL;                 /* array of parsetrees we are generating         */
   ESL_SQ      **sq_dummy = NULL;             /* dummy sequence array for output MSA           */
   CM_MX        *mx    = cm_mx_Create(cm->M); /* DP matrix for running inside algorithm        */
   int           status;                      /* esl return code                               */

   /* allocate space for parsetree array */
   ESL_ALLOC(pstr, nrow * sizeof(Parsetree_t * ));

   /* allocate space for sequence array */
   ESL_ALLOC(sq_dummy, sizeof(P7_TRACE *) * nrow);
   for (idx = 0; idx < nrow; idx++) sq_dummy[idx] = esl_sq_CreateDigital(cm->abc);

   for (n = 0; n < nseq; n++) {

      /* run inside algorithm */
      cm_InsideAlign(cm, errbuf, sq[n]->dsq, sq[n]->n, 512.0, mx, &insc);

      for (r = 0; r < R; r++) {
         idx = (n*R) + r;

         /* copy sequence into dummy array, make unique name */
         esl_sq_Copy(sq[n], sq_dummy[idx]);
         sprintf(sq_dummy[idx]->name, "%s_%d", sq[n]->name, r);

         /* sample parsetree */
         if((status = cm_StochasticParsetree(cm, errbuf, sq[n]->dsq, sq[n]->L, mx, rng, &pstr[idx], &sc)) != eslOK) return status;
         //ParsetreeDump(stdout, pstr[n*nseq + r], cm, sq[n]->dsq);
      }
   }

   /* generate MSA from sampled parsetrees */
   if ((status = Parsetrees2Alignment(cm, errbuf, cm->abc, sq_dummy, NULL, pstr, NULL, nrow, NULL, NULL, TRUE, FALSE, &msa))) return status;
   *ret_msa = msa;

   /* clean up and return */
   for (idx = 0; idx < nrow; idx++)
   {
      FreeParsetree(pstr[idx]);
      esl_sq_Destroy(sq_dummy[idx]);
   }

   free(pstr);
   cm_mx_Destroy(mx);
   free(sq_dummy);
   return eslOK;

   ERROR:
      return status;
}
