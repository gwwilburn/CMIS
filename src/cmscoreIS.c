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

#include "cmis_scoreset.h"
#include "cmis_trace.h"

/* declaration of internal functions */
int Calculate_IS_scores(CM_t *cm, P7_HMM *hmm, ESL_SQ **sq, ESL_RANDOMNESS *rng, CMIS_SCORESET *cmis_ss, ESL_MSA **msa, char *scoreprogfile, char *aliprogfile, int start, int end, int totseq, int user_R, int A, int scoreprog, int aliprog, int verbose);

static ESL_OPTIONS options[] = {
        /* name              type        default    env range togs  reqs  incomp            help                                                     docgroup */
        { "-h",              eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL,            "help; show brief info on version and usage",               1 },

        /* Options forcing which alphabet we're working in (normally autodetected) */
        { "--amino",         eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, "--dna,--rna",   "We are dealing with protein sequences",                    2 },
        { "--dna",           eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, "--amino,--rna", "We are dealing with dna sequences",                        2 },
        { "--rna",           eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, "--amino,--dna", "We are dealing with rna sequences",                        2 },

        /* options for bounding sequences we score */
        { "--seqstart",      eslARG_INT,     "-1",  NULL, NULL, NULL, NULL, NULL,            "Start sequence index",                                     3 },
        { "--seqend",        eslARG_INT,     "-2",  NULL, NULL, NULL, NULL, NULL,            "End sequence index",                                       3 },

        /* options for controlling number of alignments sampled per sequence */
        { "-R",              eslARG_INT,     "-1",  NULL, NULL, NULL, NULL, NULL,            "Number of HMM paths sampled per sequence",                 4 },
        { "-s",              eslARG_INT,     "0",   NULL, NULL, NULL, NULL, NULL,            "Set random number seed to <n>",                            4 },

        /* control of output */
        { "-A",              eslARG_OUTFILE, NULL,  NULL, NULL, NULL, NULL, NULL,            "save multiple alignment of all seqs to file <s>",          5 },
        { "--scoreprog",     eslARG_OUTFILE, NULL,  NULL, NULL, NULL, NULL, NULL,            "save .csv with info on IS scoring progress to file <s>",   5 },
        { "--aliprog",       eslARG_OUTFILE, NULL,  NULL, NULL, NULL, "-A", NULL,            "save .csv with info on IS alignment prgoress to file <s>", 5 },

        /* debugging tools */
        { "-v",              eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL,            "Verbose mode: print info on intermediate scoring steps",   6 },

        { 0,0,0,0,0,0,0,0,0,0 },
};

static char banner[] = "Align and score sequences with a CM by importance sampling a CM/SCFG";
static char usage[]  = "[-options] <cmfile> <hmmfile> <msafile> <score_outfile>";

static void
cmdline_failure(char *argv0, char *format, ...)
{
        va_list argp;
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
   char             *hmmfile          = NULL;                 /* input hmm filepath                       */
   char             *seqfile          = NULL;                 /* input seq filepath                       */
   char             *scorefile        = NULL;                 /* output score filepath                    */
   char             *scoreprogfile    = NULL;                 /* output IS scoring progress filepath      */
   char             *aliprogfile      = NULL;                 /* output IS alignment progress filepath    */
   CM_FILE          *cmfp             = NULL;                 /* open input CM file stream                */
   P7_HMMFILE       *hmmfp            = NULL;                 /* open input hmm file stream               */
   ESL_SQFILE       *sqfp             = NULL;                 /* open seq file stream                     */
   CM_t             *cm               = NULL;                 /* cm                                       */
   P7_HMM           *hmm              = NULL;                 /* hmm                                      */
   ESL_SQ          **sq               = NULL;                 /* array of sequences                       */
   CMIS_SCORESET    *cmis_ss          = NULL;                 /* for sequence scores                      */
   int               format           = eslSQFILE_UNKNOWN;    /* seq file format                          */
   FILE             *afp              = NULL;                 /* output alignment file (-A)               */
   ESL_MSA          *msa              = NULL;                 /* alignment output object (-A)             */
   int               outfmt           = eslMSAFILE_STOCKHOLM; /* alignment output format (-A)             */
   FILE             *cmis_ss_fp       = NULL;                 /* file for sequence scores                 */
   int               n                =  0;                   /* seq_index                                */
   int               nseq             =  0;                   /* number of seqs in seq file               */
   int               totseq           =  0;                   /* number of seqs we deal with              */
   int               start,end;                               /* start and end seq indices                */
   int               R                = -1;                   /* number of ali samples/seq                */
   int               A                = 0;                    /* Boolean for output alignment             */
   int               scoreprog        = 0;                    /* Boolean for IS scoring progress output   */
   int               aliprog          = 0;                    /* Boolean for IS alignment progress output */
   int               v                = 0;                    /* Boolean for verbose mode                 */
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
      puts("\n Options for bounding sequences we score:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
      puts("\n Sampling options:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80);
      puts("\n MSA Output options:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80);
      puts("\n Debug options:");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80);
      exit(0);
   }

   /* read arguments */
   if (esl_opt_ArgNumber(go) != 4)  cmdline_failure(argv[0], "Incorrect number of command line arguments.\n", go->errbuf);

   cmfile    =  esl_opt_GetArg(go, 1);
   hmmfile   =  esl_opt_GetArg(go, 2);
   seqfile   =  esl_opt_GetArg(go, 3);
   scorefile =  esl_opt_GetArg(go, 4);

   /* check for verbose mode */
   if (esl_opt_GetBoolean(go, "-v")) v=1;

   /* if output msa requested by user, try to open it */
   if (esl_opt_IsOn(go, "-A")) {
      if ((afp = fopen(esl_opt_GetString(go, "-A"), "w")) == NULL) p7_Fail("Failed to open alignment file %s for writing\n", esl_opt_GetString(go, "-A"));
      A = 1;
   }

   /* if IS intermediate scoring progress csv is requested by user, set variables*/
   if (esl_opt_IsOn(go, "--scoreprog")) {
      scoreprogfile = (esl_opt_GetString(go, "--scoreprog"));
      scoreprog = 1;
   }

   /* if IS intermediate alignment progress csv is requested by user, set variables*/
   if (esl_opt_IsOn(go, "--aliprog")) {
      aliprogfile = (esl_opt_GetString(go, "--aliprog"));
      aliprog = 1;
   }


   /* if user has defined an alphabet we define it here */
   if        (esl_opt_GetBoolean(go, "--amino"))       abc = esl_alphabet_Create(eslAMINO);
   else if   (esl_opt_GetBoolean(go, "--rna"))         abc = esl_alphabet_Create(eslRNA);
   else if   (esl_opt_GetBoolean(go, "--dna"))         abc = esl_alphabet_Create(eslDNA);

   /* check if user has specified the number of samples/sequence */
   if (esl_opt_GetInteger(go, "-R") > 0) R = esl_opt_GetInteger(go, "-R");

   /* create random number generator */
   rng = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

   /* open the .hmm file */
   status = p7_hmmfile_OpenE(hmmfile, NULL, &hmmfp, errbuf);
   if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", hmmfile, errbuf);
   else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                hmmfile, errbuf);
   else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, hmmfile, errbuf);

   /* read first hmm  from hmm file */
   if (p7_hmmfile_Read(hmmfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

   p7_hmmfile_Close(hmmfp);

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
   cm_Configure(cm, errbuf, -1);

   /* open the sequence file */
   status = esl_sqfile_OpenDigital(abc, seqfile, format, NULL, &sqfp);
   if      (status == eslENOTFOUND) p7_Fail("No such file.");
   else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
   else if (status ==eslEINVAL)     p7_Fail("Can't autodetect stdin or .gz.");
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

   /* if user specified bounds, set them here */
   start = 0;
   end   = nseq;
   if (esl_opt_GetInteger(go, "--seqstart") > -1) start = esl_opt_GetInteger(go, "--seqstart");
   if (esl_opt_GetInteger(go, "--seqend")   > -1)     end = esl_opt_GetInteger(go, "--seqend");
   /* check if end is beyond number of sequences in input file */
   if (end > nseq) end = nseq;


   /* total number of sequences we will be scoring */
   totseq = end-start;
   fprintf(stdout, "totseq: %d\n", totseq);

   /* create scoreset object */
   cmis_ss = cmis_scoreset_Create(totseq);

   /* calculate importance sampling score for all seqs */
   Calculate_IS_scores(cm, hmm, sq, rng, cmis_ss, &msa, scoreprogfile, aliprogfile, start, end, totseq, R, A, scoreprog, aliprog, v);

   /* write score file */
   if ((cmis_ss_fp = fopen(scorefile, "w")) == NULL) esl_fatal("Failed to open output CMIS scoreset file %s for writing", scorefile);
   if (cmis_scoreset_Write(cmis_ss_fp, cmis_ss) != eslOK) esl_fatal("Failed to write scores to CMIS scoreset file %s",    scorefile);
   fclose(cmis_ss_fp);

   if(cm->flags & CMH_LOCAL_BEGIN) fprintf(stdout, "Local begins on!\n");
   if(cm->flags & CMH_LOCAL_END) fprintf(stdout, "Local ends on!\n");

   /* if output msa requested, write it to output file */
   if (A) {
      esl_msafile_Write(afp, msa, outfmt);
      fclose(afp);
      esl_msa_Destroy(msa);
   }

   /* clean up and return */
   for (n = 0; n < nseq+1; n++)
   {
      esl_sq_Destroy(sq[n]);
   }

   esl_sqfile_Close(sqfp);
   free(sq);
   FreeCM(cm);
   cmis_scoreset_Destroy(cmis_ss);
   p7_hmm_Destroy(hmm);
   esl_randomness_Destroy(rng);
   esl_alphabet_Destroy(abc);
   esl_getopts_Destroy(go);
   return 0;

   ERROR:
   return status;
}


int Calculate_IS_scores(CM_t *cm, P7_HMM *hmm, ESL_SQ **sq, ESL_RANDOMNESS *rng, CMIS_SCORESET *cmis_ss, ESL_MSA **msa, char *scoreprogfile, char *aliprogfile, int start, int end, int totseq, int user_R, int A, int scoreprog, int aliprog, int verbose)
{
   P7_PROFILE     *gm          = NULL;                      /* h3 profile model                            */
   P7_BG          *bg          = NULL;                      /* h3 background model                         */
   P7_GMX         *fwd         = p7_gmx_Create(100, 100);   /* h3 forward matrix                           */
   P7_TRACE      **tr_dummy;                                /* dummy trace array for sampled alignments    */
   ESL_SQ        **sq_dummy;                                /* dummy sequence array for samp;ed alignments */
   int             msaopts     = 0;                         /* options for creating an MSA                 */
   P7_TRACE      **out_tr      = NULL;                      /* traces used to construct output msa (-A)    */
   ESL_SQ        **out_sq      = NULL;                      /* traces used to construct output msa (-A)    */
   FILE           *scoreprogfp = NULL;                      /* score progress output file (--scoreprog)    */
   FILE           *aliprogfp   = NULL;                      /* alignment progress output file (--aliprog)  */
   int             i,j,k;                                   /* sequence indices                            */
   int             r;                                       /* alignment sample index                      */
   int             cpos;                                    /* cm node index                               */
   int             R           = 1000;                      /* number of sampled alignments per sequence   */
   int             nBetter;                                 /* number of better alignments we've found     */
   float           fsc;                                     /* forward partial log-odds score, in nats     */
   float           ntsc;                                    /* hmmer null transition score, in nats        */
   float           hmmsc_ld;                                /* hmm log odds S(x,pi), in nats               */
   float           cmsc_ld;                                 /* cm log odds S(x, pi), in bits               */
   float           cmsc_max;                                /* best cm log odds S(x, pi) for all paths     */
   float           insc;                                    /* log-odds score for inside algorithm         */
   float           ls;                                      /* log of importance sampling sum              */
   float           ld;                                      /* importance sampling log odds score S(x)     */
   int             status;                                  /* easel return code                           */
   char            errbuf[eslERRBUFSIZE];                   /* buffer for easel errors                     */


   /* if score progress file is requested, open it */
   if (scoreprog) {
      if ((scoreprogfp = fopen(scoreprogfile, "w")) == NULL) p7_Fail("Failed to open file %s for writing\n", scoreprogfile);
      fprintf(scoreprogfp, "id,iter,inside_logodds,cmis_logodds\n");
   }

   /* if alignment progress file is requested, open it */
   if (aliprog) {
      if ((aliprogfp = fopen(aliprogfile, "w")) == NULL) p7_Fail("Failed to open file %s for writing\n", aliprogfile);
      fprintf(aliprogfp, "id,iter,cmis_logodds\n");
   }

   /* set number of samples */
   if (user_R > 0) R = user_R;

   /* allocate space for dummy trace  and sequence arrays */
   ESL_ALLOC(tr_dummy, sizeof(P7_TRACE *) * R);
   ESL_ALLOC(sq_dummy, sizeof(P7_TRACE *) * R+1);
   for (r = 0; r < R; r++) {
      tr_dummy[r] = p7_trace_Create();
      sq_dummy[r] = esl_sq_CreateDigital(hmm->abc);
   }

   if (A) {
      /* allocate memory for traces  and initialize */
      ESL_ALLOC(out_tr, sizeof(P7_TRACE *) * totseq);
      ESL_ALLOC(out_sq, sizeof(ESL_SQ *) * (totseq+1));
   }

   /* create a profile from the HMM */
   gm = p7_profile_Create(hmm->M, hmm->abc);

   /* create null model*/
   bg = p7_bg_Create(hmm->abc);

   /* configure profile with dummy length (we'll change it in outer for loop) */
   p7_ProfileConfig(hmm, bg, gm, sq[0]->n, p7_UNIGLOCAL);

   /* set up cm for scoring */
   if((status = CMLogoddsify(cm)) != eslOK) ESL_FAIL(status, errbuf, "problem logodisfying CM");


   /* outer loop over sequences */
   for (i=start; i < end; i++) {

      nBetter = 0;

      /* index for score set: must start at 0 */
      j=i-start;

      /* for getting optimal path under hpm */
      if (A) {
         cmsc_max  = -eslINFINITY;
         out_sq[j] = esl_sq_CreateDigital(hmm->abc);
         esl_sq_Copy(sq[i], out_sq[j]);
      }



      /* declare variables that are re-initialized at each iteration of outer loop */
      ESL_MSA        *msa_dummy;                            /* MSA of different aligmnetss of seq i  */
      Parsetree_t    *mtr       = NULL;                     /* the guide tree                        */
      Parsetree_t   **pstr      = NULL;                     /* array of parse trees created from MSA */
      CM_MX          *mx        = cm_mx_Create(cm->M);      /* DP matrix for inside algorithm)       */

      /* Set the profile and null model's target length models */
      p7_ReconfigLength(gm, sq[i]->n);
      p7_bg_SetLength(bg, sq[i]->n);

      /* resize dp matrix if necessary */
      p7_gmx_GrowTo(fwd, gm->M, sq[i]->n);

      /* run forward algorithm */
      p7_GForward(sq[i]->dsq, sq[i]->n, gm, fwd, &fsc);

      /* calculate null transition score */
      p7_bg_NullOne(bg, sq[i]->dsq, sq[i]->n, &ntsc);

      /* calculate Inside log odds score */
      cm_InsideAlign(cm, errbuf, sq[i]->dsq, sq[i]->n, 512.0, mx, &insc);

      /* allocate variables to store sequence-alignment log-odds scores */
      double *pr;
      ESL_ALLOC(pr, R*sizeof(double));
      esl_vec_DSet(pr, R, 0.0);



      /* inner loop 1: sample alignments from hmm */
      for (r = 0; r < R; r++) {

         /* copy sequence into dummy array */
         esl_sq_Copy(sq[i], sq_dummy[r]);

         /* get stochastic traceback */
         p7_GStochasticTrace(rng, sq[i]->dsq, sq[i]->n, gm, fwd, tr_dummy[r]);

         //p7_trace_Dump(stdout, tr_dummy[r], gm, sq[i]->dsq);

         /* get HMM score S(x,pi) for this seq, trace */
         p7_trace_Score(tr_dummy[r], sq[i]->dsq, gm, &hmmsc_ld);

         /* convert score to bits, put in  array */
         pr[r] -= (hmmsc_ld) /  eslCONST_LOG2;

      }

      /* convert traces to ESL_MSA */
      msaopts |= p7_ALL_CONSENSUS_COLS; /* keep all consensus columns                 */
      msaopts |= p7_DIGITIZE;           /* make digital msa (required for guide tree) */
      p7_tracealign_Seqs(sq_dummy, tr_dummy, R, hmm->M, msaopts, hmm, &msa_dummy);

      /* map cm's SS_cons to MSA */
      ESL_ALLOC(msa_dummy->ss_cons, (sizeof(char) * (msa_dummy->alen+1)));
      cpos = 0;
      for (k = 0; k <= msa_dummy->alen; k++) {
         if (msa_dummy->rf[k] == 120) {
            msa_dummy->ss_cons[k] = cm->cmcons->cstr[cpos];
            cpos++;
         }
         else {
            msa_dummy->ss_cons[k] = '.';
         }
      }
      msa_dummy->ss_cons[msa_dummy->alen] = '\0';

      /* create guide tree */
      status = HandModelmaker(msa_dummy, errbuf,
                              TRUE,  /* use_rf */
                              FALSE, /* use_el, no */
                              FALSE, /* use_wts, irrelevant */
                              0.5,   /* gapthresh, irrelevant */
                              NULL,  /* returned CM, irrelevant */
                              &mtr); /* guide tree */
      if (status != eslOK) return status;

      /* create parse tree alignment from MSA */
      Alignment2Parsetrees(msa_dummy, cm, mtr, errbuf, NULL, &pstr);

      /* inner loop 2: score traces with CM */
      for (r = 0; r < R; r++) {

         //ParsetreeDump(stdout, pstr[r], cm, sq[i]->dsq);

         /* score this sequence and alignment */
         ParsetreeScore(cm,
                        NULL,
                        errbuf,
                        pstr[r],
                        sq[i]->dsq,
                        0,
                        &cmsc_ld,
                        NULL,
                        NULL,
                        NULL,
                        NULL);

         pr[r] += cmsc_ld;

         /* if requested, track intermediate scoring progress */
         if (( (r+1)  % 10000 == 0) && (scoreprog)) {

            esl_vec_DSortIncreasing(pr, r+1);
            float ls = esl_vec_DLog2Sum(pr, r+1);
            float ld = ((fsc - logf(r+1) ) / eslCONST_LOG2) + ls;


            //fprintf(stdout, "\n\n\t## intermediate details for sequence %s, iteration %d\n", sq[i]->name, r);
            fprintf(scoreprogfp, "%s,%d,%.2f,%.2f\n", sq[i]->name, r, insc, ld);


         }

         /* for output alignment , see if we've got a better path */
         if (A && cmsc_ld > cmsc_max) {
            cmsc_max = cmsc_ld;

            /* if we're tracking alignment progress, print S(x,pi) of new best alignment */
            if (aliprog) fprintf(aliprogfp, "%s,%d,%.2f\n", sq[i]->name, r, cmsc_max);

            /* if we've already created out_tr[j], free it before recreating it */
            if (nBetter > 0 ) p7_trace_Destroy(out_tr[j]);
            out_tr[j] = cmis_trace_Clone(tr_dummy[r]);
            nBetter++;
         }
         //fprintf(stdout, "\talignment: %d, hmmsc_ld: %.4f, cmsc_ld: %.4f\n,", r, hmmsc_ld /  eslCONST_LOG2, cmsc_ld);

         //fprintf(stdout, "CM score for sequence %d and alignment %d: %.4f\n", i, r, cmsc_ld);
         /* clean up for next iteration of outer loop */
         p7_trace_Reuse(tr_dummy[r]);
         esl_sq_Reuse(sq_dummy[r]);
         FreeParsetree(pstr[r]);
      }


      /* calculate IS approximation for this sequence */

      /* sort to increase accuracy if LogSum */
      esl_vec_DSortIncreasing(pr, R);

      /* run log sum */
      ls = esl_vec_DLog2Sum(pr, R);

      /* incorporate alignment-independent bits into score */
      ld = ((fsc - logf(R) ) / eslCONST_LOG2) + ls;

      /* add scoring info to scoreset object */
      cmis_ss->sqname[j]  = sq[i]->name;
      cmis_ss->R[j]       = R;
      cmis_ss->fsc[j]     = fsc / eslCONST_LOG2;
      cmis_ss->ntsc[j]    = ntsc / eslCONST_LOG2;
      cmis_ss->insc[j]    = insc;
      cmis_ss->cmis_ld[j] = ld;


      /* clean up for next iteration of outer loop */

      esl_msa_Destroy(msa_dummy);
      free(pstr);
      FreeParsetree(mtr);
      cm_mx_Destroy(mx);
      free(pr);

   }

   /* if output msa requested, create it from traces */
   if (A) {
      p7_tracealign_Seqs(out_sq, out_tr, totseq, hmm->M, msaopts, hmm, msa);
   }

   /* clean up and return */
   if (scoreprogfp) fclose(scoreprogfp);
   if (aliprogfp) fclose(aliprogfp);
   if (A) {
      for (j=0; j<totseq; j++) {
         p7_trace_Destroy(out_tr[j]);
         esl_sq_Destroy(out_sq[j]);
      }
      free(out_tr);
      free(out_sq);
   }

   for (r = 0; r < R; r++) {
      esl_sq_Destroy(sq_dummy[r]);
   }
   free(sq_dummy);
   p7_trace_DestroyArray(tr_dummy, R);
   p7_profile_Destroy(gm);
   p7_bg_Destroy(bg);
   p7_gmx_Destroy(fwd);
   return eslOK;

   ERROR:
      return status;
}

