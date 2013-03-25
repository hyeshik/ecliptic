/*
 * crossfest.c
 *   - Generates simulated random profiles based on experimental parameters
 *
 * Copyright (c) 2012 Hyeshik Chang.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <zlib.h>
#include <sys/param.h>
#include <sys/time.h>
#include <pthread.h>
#include "clipstatstree.h"

#define SFMT_MEXP       19937
#include "SFMT.h"

#define TIMEVAL_DIFF(b, a) \
    ((double)((a).tv_sec - (b).tv_sec) + ((a).tv_usec - (b).tv_usec) * 0.000001)

#define MAXSEQLEN       8388608 /* change this if you see a larger cluster */
#define MAXFRAGNUM      32      /* change this if more spliced reads are read */

#define NUMBASES        5
#define BASE_A          0
#define BASE_C          1
#define BASE_G          2
#define BASE_T          3
#define BASEDELETION    4

static const int nucleobase2int[256] = { 
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 4, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 0, 9, 1, 9, 9, 9, 2, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 3, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 0, 9, 1, 9, 9, 9, 2, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 3, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
};
static const char *NUCLEOBASES="ACGT-";

typedef struct {
    uint32_t readlength;
    uint32_t numbasetypes;
} ERROR_PROFILE_HEADER;

typedef struct {
    uint32_t readlength;
    uint32_t _padding;
    uint64_t readdist[1][NUMBASES][NUMBASES];
} ERROR_PROFILE;

typedef struct {
    unsigned char *reads;
    unsigned char *sentinel;
    unsigned char *next;
    void *padding; /* padding */
} READ_QUEUE;

typedef struct {
    uint32_t readlength;
    READ_QUEUE queues[1][NUMBASES];
} READ_QUEUE_SET;
#define READ_QUEUE_UNIFORM_EXTENSION    10  /* for prolonged reads including */
                                            /* deletion */

typedef struct {
    char name[32];          /* reference sequence name (eg. chr1) */
    uint32_t start;         /* start coordinate of this window (- for neg strand) */
    uint32_t end;           /* end coordinate of this window (- for neg strand) */
    uint32_t nuniquereads;  /* number of unique reads (# of following records) */
    uint32_t seqlength;     /* sequence length in bytes of this window (=nucleotides) */
} READCLUSTER_HEADER;

#pragma pack(1)
typedef struct {
    uint32_t nreads;        /* number of duplicate reads (instances) */
    uint8_t nfragments;     /* number of fragments (>= 2 for spliced exons) */
} MAPPING_HEADER;
#pragma pack()

typedef struct {
    uint32_t start;
    uint32_t end;
} MAPPING_FRAGMENT;

typedef struct {
    uint32_t total;
    uint32_t done;
    uint32_t queued;
    int nthreads;
    struct timeval started;
    pthread_mutex_t lock;
} JOB_COUNTER;

struct _worker;
typedef struct _worker {
    struct _worker *workers;
    pthread_t thread;
    ERROR_PROFILE *error_profile;
    const char *readpool_input_path;
    JOB_COUNTER *jobcounter;
    int threadid;
    double last_consumed;
    struct timeval last_started;
    int status;
    int show_progression;
} WORKER;
#define WORKER_STATUS_NOT_RUNNING           '/'
#define WORKER_STATUS_READ_QUEUE_SHUFFLING  'r'
#define WORKER_STATUS_SIMULATING            's'
#define WORKER_STATUS_MERGING_OUTPUT        'M'
#define WORKER_STATUS_FINISHED              '_'
#define WORKER_STATUS_ERROR                 'X'
#define WORKER_STATUS_IS_RUNNING(s)                 \
    ((s) == WORKER_STATUS_READ_QUEUE_SHUFFLING ||   \
     (s) == WORKER_STATUS_SIMULATING)

#define RANDGEN_BLOCKSIZE   16384
/* NOT THREAD-SAFE IN FAVOR OF SPEED. ALLOCATE THIS PER-THREAD BASIS */
typedef struct {
    sfmt_t sfmt;
    uint64_t rarray[RANDGEN_BLOCKSIZE];
    int rarraynext;
} RANDGEN_STATE;

void
usage(const char *execpath)
{
    printf("CROSSFEST version 1.1 - crosslinking FDR estimator\n");
    printf("Part of `nrclip' by Hyeshik Chang <hyeshik@snu.ac.kr>\n\n");
    printf("This program implements a permutation-based approach described "
           "in Zhang et al. (Nature Biotechnology, 2011, 29(7):607-14).\n\n");
    printf("Usage: %s [OPTIONS]\n\n", execpath);
    printf("Mapping and error profile files are prepared by "
           "an external script.\n\n");
    printf("  -i FILE\tinput readpool file (required)\n");
    printf("  -e FILE\tinput error profile array (required)\n");
    printf("  -o FILE\toutput file prefix (required)\n");
    printf("  -t N\t\tnumber of worker threads (default: 8)\n");
    printf("  -d N\t\tminimum depth of simulated reads in Gb (default: 500)\n");
}

static void
randgen_init(RANDGEN_STATE *rstate)
{
    sfmt_init_gen_rand(&rstate->sfmt, rand());
    sfmt_fill_array64(&rstate->sfmt, rstate->rarray, RANDGEN_BLOCKSIZE);
    rstate->rarraynext = 0;
}

/* `maximum' must be far smaller than 2^32 to minimize modulo biases */
static inline uint64_t
randgen_rand(RANDGEN_STATE *rstate, uint64_t maximum)
{
    if (rstate->rarraynext >= RANDGEN_BLOCKSIZE) {
        sfmt_fill_array64(&rstate->sfmt, rstate->rarray, RANDGEN_BLOCKSIZE);
        rstate->rarraynext = 0;
    }

    return rstate->rarray[rstate->rarraynext++] % (maximum + 1);
}

static ERROR_PROFILE *
load_error_profile(const char *filepath)
{
    FILE *fp=NULL;
    ERROR_PROFILE_HEADER filehead;
    ERROR_PROFILE *profile=NULL;
    size_t numrecords, recsize;

    fp = fopen(filepath, "r");
    if (fp == NULL) {
        fprintf(stderr, "Failed to open file %s\n", filepath);
        goto onError;
    }

    if (fread(&filehead, sizeof(filehead), 1, fp) < 1) {
        fprintf(stderr, "Unexpected end of error profile.\n");
        goto onError;
    }

    if (filehead.numbasetypes != NUMBASES) {
        fprintf(stderr, "Number of base types is %d, which is not supported.\n",
                filehead.numbasetypes);
        goto onError;
    }

    numrecords = filehead.readlength * filehead.numbasetypes;
    recsize = sizeof(uint64_t) * NUMBASES;

    profile = (ERROR_PROFILE *)malloc(sizeof(uint32_t) * 2 +
                                      recsize * numrecords);
    if (profile == NULL) {
        fprintf(stderr, "Can't allocate memory for error profile.\n");
        goto onError;
    }

    if (fread(profile->readdist, recsize, numrecords, fp) != numrecords) {
        fprintf(stderr, "Incomplete error profile.\n");
        goto onError;
    }

    profile->readlength = filehead.readlength;
    printf("Loaded error profile for %d cycles.\n", filehead.readlength);

    fclose(fp);

    return profile;

  onError:
    if (fp != NULL)
        fclose(fp);
    if (profile != NULL)
        free(profile);

    return NULL;
}

static int
generate_shuffled_read_queue(RANDGEN_STATE *rgen, READ_QUEUE *queue,
                             const uint64_t *targetbases,
                             unsigned char defaultmark)
{
    uint64_t totalreads=0, i;
    unsigned char *cur;
    int actualbase;

    for (actualbase = 0; actualbase < NUMBASES; actualbase++)
        totalreads += targetbases[actualbase];

    /* Allocate space and prepare a repetitive string represents frequency of
     * actual read bases. */
    queue->reads = queue->next = malloc(totalreads > 1 ? totalreads : 1);
    if (queue->reads == NULL)
        return -1;

    if (totalreads > 0)
        queue->sentinel = queue->reads + totalreads;
    else {
        /* Set a default base for blank queues */
        queue->reads[0] = defaultmark;
        queue->sentinel = queue->reads + 1;
        return 0;
    }

    cur = queue->reads;
    for (actualbase = 0; actualbase < NUMBASES; actualbase++) {
        memset(cur, actualbase, targetbases[actualbase]);
        cur += targetbases[actualbase];
    }

    assert(cur == queue->sentinel);

    /* Fisher-Yates shuffle of the read bases */
    if (totalreads > 1) {
        for (i = totalreads - 1; i >= 1; i--) {
            uint64_t j;
            unsigned char temp;

            j = randgen_rand(rgen, i);
            temp = queue->reads[i];
            queue->reads[i] = queue->reads[j];
            queue->reads[j] = temp;
        }
    }

    return 0;
}

static void
free_read_queue_set(READ_QUEUE_SET *queueset)
{
    int pos, refbase;

    for (pos = 0; pos < queueset->readlength; pos++)
        for (refbase = 0; refbase < NUMBASES; refbase++)
            if (queueset->queues[pos][refbase].reads != NULL)
                free(queueset->queues[pos][refbase].reads);

    free(queueset);
}

static READ_QUEUE_SET *
generate_read_queue_set(RANDGEN_STATE *rgen, ERROR_PROFILE *error_profile)
{
    READ_QUEUE_SET *queueset;
    size_t numqueues;
    uint32_t pos, refbase, run_readlength;
    const uint64_t dummy_readdist[NUMBASES]={0, 0, 0, 0, 0};

    numqueues = (error_profile->readlength + READ_QUEUE_UNIFORM_EXTENSION)
                * NUMBASES;
    queueset = (READ_QUEUE_SET *)malloc(sizeof(READ_QUEUE) * numqueues);
    if (queueset == NULL) {
        fprintf(stderr, "Failed to allocate memory for random read queue.\n");
        return NULL;
    }

    memset(queueset, 0, sizeof(READ_QUEUE) * numqueues);
    run_readlength = error_profile->readlength;

    for (pos = 0; pos < run_readlength; pos++) {
        for (refbase = 0; refbase < NUMBASES; refbase++) {
            if (generate_shuffled_read_queue(rgen,
                    &queueset->queues[pos][refbase],
                    error_profile->readdist[pos][refbase], refbase) != 0) {
                fprintf(stderr, "Error on shuffling read queue\n");
                goto onError;
            }
        }
    }

    /* Extended uniform read queues for prolonged reads by deletion */
    for (pos = 0; pos < READ_QUEUE_UNIFORM_EXTENSION; pos++) {
        for (refbase = 0; refbase < NUMBASES; refbase++) {
            if (generate_shuffled_read_queue(rgen,
                    &queueset->queues[pos + run_readlength][refbase],
                    dummy_readdist, refbase) != 0) {
                fprintf(stderr, "Error on padding extended read queue\n");
                goto onError;
            }
        }
    }
    queueset->readlength = run_readlength + READ_QUEUE_UNIFORM_EXTENSION;

    return queueset;

  onError:
    free_read_queue_set(queueset);
    return NULL;
}

static void
convert_nucleotide_seq_to_num(unsigned char *seq, size_t size)
{
    unsigned char *cur;

    for (cur = seq + size - 1; cur >= seq; cur--)
        *cur = nucleobase2int[(int)*cur];
}

static inline int
check_if_sequence_valid(unsigned char *seq, size_t start, size_t end)
{
    size_t pos;

    for (pos = start; pos < end; pos++)
        if (seq[pos] >= NUMBASES)
            return 0;

    return 1;
}

static double
shannon_entropy(const uint32_t *counts)
{
    ssize_t i;
    uint32_t maxreads, totalreads;
    double entropy;
    double totalreadsf;
    maxreads = totalreads = 0;

    for (i = 0; i < NUMBASES; i++) {
        totalreads += counts[i];
        if (maxreads < counts[i])
            maxreads = counts[i];
    }   

    if (maxreads == totalreads)
        return 0.; 

    totalreadsf = (double)totalreads;
    entropy = 0.; 
    for (i = 0; i < NUMBASES; i++) {
        double p = counts[i] / totalreadsf;
        if (p > 0.)
            entropy -= log(p) * p;
    }

    return entropy;
}

static int
simulate_sequencing(gzFile *inputf, READ_QUEUE_SET *queueset, CLIPSTATS_TREES *trees)
{
    READCLUSTER_HEADER header;
    MAPPING_FRAGMENT *mapping=NULL;
    unsigned char *refseq=NULL;
    uint32_t *posreadcount=NULL;

    refseq = malloc(MAXSEQLEN);
    if (refseq == NULL) {
        fprintf(stderr, "Unable to allocate sequence buffer\n");
        goto onError;
    }

    posreadcount = malloc(sizeof(uint32_t) * MAXSEQLEN * NUMBASES);
    if (posreadcount == NULL) {
        fprintf(stderr, "Unable to allocate read count buffer\n");
        goto onError;
    }

    mapping = malloc(sizeof(MAPPING_FRAGMENT) * MAXFRAGNUM);
    if (mapping == NULL) {
        fprintf(stderr, "Unable to allocate mapping buffer\n");
        goto onError;
    }

    for (;;) {
        int r;
        uint32_t i;
        uint32_t *cntptr;

        r = gzread(inputf, &header, sizeof(header));
        if (r == -1) {
            fprintf(stderr, "Error occurred while reading mapping file.\n");
            goto onError;
        }
        else if (r == 0)
            break;

        if (header.seqlength > MAXSEQLEN) {
            fprintf(stderr, "Sequence `%s:%u-%u' too long. Increase MAXSEQLEN in the code "
                    "and compile, re-run the analysis.\n",
                    header.name, header.start, header.end);
            goto onError;
        }

        if (gzread(inputf, refseq, header.seqlength) < header.seqlength) {
            fprintf(stderr, "Failed to load sequence `%s:%u-%u'\n", header.name,
                    header.start, header.end);
            goto onError;
        }

        convert_nucleotide_seq_to_num(refseq, header.seqlength);
        memset(posreadcount, 0, sizeof(uint32_t) * header.seqlength * NUMBASES);

        /* Simulate sequencing of reads for a transcript */
        for (i = 0; i < header.nuniquereads; i++) {
            MAPPING_HEADER mappingheader;
            uint32_t n;
            int j, k;

            if (gzread(inputf, &mappingheader, sizeof(mappingheader)) <
                    sizeof(mappingheader)) {
                fprintf(stderr, "Failed to load a mapping in `%s:%u-%u'\n", header.name,
                        header.start, header.end);
                goto onError;
            }

            if (mappingheader.nfragments > MAXFRAGNUM) {
                fprintf(stderr, "Too many fragments (%d) for a mapping in `%s:%u-%u'\n",
                        mappingheader.nfragments, header.name, header.start, header.end);
                continue;
            }

            if (gzread(inputf, mapping, sizeof(MAPPING_FRAGMENT) * mappingheader.nfragments) <
                       sizeof(MAPPING_FRAGMENT) * mappingheader.nfragments) {
                fprintf(stderr, "Failed to load a mapping list in `%s:%u-%u'\n", header.name,
                        header.start, header.end);
                goto onError;
            }

            for (n = 0; n < mappingheader.nreads; n++) {
                int fragi;

                /* each fragment in a unique read (insertion or deletion may push or
                 * pull read positions) */
                k = 0;

                for (fragi = 0; fragi < mappingheader.nfragments; fragi++) {
                    /* j for refseq position, k for read position */
                    int fragend=mapping[fragi].end;

                    for (j = mapping[fragi].start;
                            j < fragend && k < queueset->readlength; j++, k++) {
                        int readgen, refbase=refseq[j];
                        READ_QUEUE *queue;

                        queue = &queueset->queues[k][refbase];
                        readgen = *(queue->next++);
                        if (queue->next >= queue->sentinel)
                            queue->next = queue->reads; /* rewind the queue */

                        posreadcount[j * NUMBASES + readgen]++;

                        /* lagging refseq index for deletions */
                        if (readgen == BASEDELETION)
                            k--;
                    }
                }
            }
        }

        /* Calculate reporting measures */
        cntptr = posreadcount;
        for (i = 0; i < header.seqlength; i++, cntptr += NUMBASES) {
            int j;
            uint32_t postotalreads, refbaseread;
            double delrate, modrate, moddelrate, entropy, t2crate;
            uint8_t basetype;

            postotalreads = 0;
            for (j = 0; j < NUMBASES; j++)
                postotalreads += cntptr[j];

            if (postotalreads < 1)
                continue;

            delrate = cntptr[BASEDELETION] / (double)postotalreads;
            refbaseread = cntptr[refseq[i]];
            modrate = (postotalreads - refbaseread - cntptr[BASEDELETION])
                        / (double)postotalreads;
            moddelrate = (postotalreads - refbaseread) / (double)postotalreads;
            entropy = shannon_entropy(cntptr);
            basetype = refseq[i];
            t2crate = (basetype == BASE_T) ? (cntptr[BASE_C] / (double)postotalreads) : 0.;

            clipstats_add(trees, &trees->del, postotalreads, basetype, delrate);
            clipstats_add(trees, &trees->mod, postotalreads, basetype, modrate);
            clipstats_add(trees, &trees->moddel, postotalreads, basetype, moddelrate);
            clipstats_add(trees, &trees->entropy, postotalreads, basetype, entropy);
            clipstats_add(trees, &trees->t2crate, postotalreads, basetype, t2crate);
        }
    }

    free(mapping);
    free(posreadcount);
    free(refseq);

    return 0;

  onError:
    if (mapping != NULL)
        free(mapping);
    if (refseq != NULL)
        free(refseq);
    if (posreadcount != NULL)
        free(posreadcount);

    return -1;
}

/* Assume the job counter mutex is locked inside this function */
static void
update_progression(WORKER *worker)
{
    JOB_COUNTER *jc;
    double total_job_time;
    int i, ran_threads;

    jc = worker->jobcounter;

    ran_threads = 0;
    total_job_time = 0.;

    if (worker->show_progression)
        printf("\r%5.1f%% done  [", (100. * jc->done) / jc->total);

    for (i = 0; i < jc->nthreads; i++) {
        if (worker->show_progression)
            printf("%c", worker->workers[i].status);
        if (worker->workers[i].last_consumed >= 0.) {
            total_job_time += worker->workers[i].last_consumed;
            ran_threads++;
        }
    }
    if (worker->show_progression)
        printf("]");

    if (ran_threads >= 1) {
        double mean_job_time, eta_waiting_jobs, eta_running_jobs, eta;
        struct timeval tv;

        gettimeofday(&tv, NULL);
        mean_job_time = total_job_time / ran_threads;

        eta_waiting_jobs = ceil(((double)(jc->total - jc->queued))
                                / jc->nthreads) * mean_job_time;
        eta_running_jobs = 0.;
        for (i = 0; i < jc->nthreads; i++)
            if (WORKER_STATUS_IS_RUNNING(worker->workers[i].status)) {
                double remaining;

                remaining = mean_job_time - TIMEVAL_DIFF(
                                worker->workers[i].last_started, tv);
                if (remaining < 0.)
                    remaining = 0.;
                if (remaining > eta_running_jobs)
                    eta_running_jobs = remaining;
            }

        eta = eta_running_jobs + eta_waiting_jobs;

        if (eta > 0) {
            time_t eft_timestamp;
            struct tm *eft_tm;
            char eft_asc[BUFSIZ];

            time(&eft_timestamp);
            eft_timestamp += eta;
            eft_tm = localtime(&eft_timestamp);
            strftime(eft_asc, BUFSIZ-1, "%b %d %H:%M:%S", eft_tm);
            if (worker->show_progression)
                printf("  Est.Fin. %s   ", eft_asc);
        }
    }

    if (worker->show_progression)
        fflush(stdout);
}

static void *
run_clip_permutation(void *args)
{
    CLIPSTATS_TREES *trees;
    JOB_COUNTER *jobcounter;
    RANDGEN_STATE randgen;
    WORKER *worker;
    gzFile fp=NULL;

    randgen_init(&randgen);
    worker = (WORKER *)args;
    jobcounter = worker->jobcounter;

    fp = gzopen(worker->readpool_input_path, "r"); 
    if (fp == NULL) {
        fprintf(stderr, "Failed to open input file %s\n", worker->readpool_input_path);
        return NULL;
    }

    trees = clipstatstrees_new();

#define BEGIN_JOBCOUNTER_EXCLUSIVE  pthread_mutex_lock(&jobcounter->lock);
#define END_JOBCOUNTER_EXCLUSIVE    pthread_mutex_unlock(&jobcounter->lock);
#define LOCK_AND_UPDATE_STATUS(worker, newstatus) do {  \
    pthread_mutex_lock(&(worker)->jobcounter->lock);    \
    (worker)->status = newstatus;                       \
    update_progression((worker));                       \
    pthread_mutex_unlock(&(worker)->jobcounter->lock);  \
} while (0)

    for (;;) {
        READ_QUEUE_SET *read_queue_set=NULL;
        struct timeval tv;

        /* Increase the job counter, exit the loop if all works done */
        BEGIN_JOBCOUNTER_EXCLUSIVE

            if (jobcounter->queued >= jobcounter->total) {
                END_JOBCOUNTER_EXCLUSIVE
                break;
            }

            jobcounter->queued++;
            gettimeofday(&worker->last_started, NULL);
            worker->status = WORKER_STATUS_READ_QUEUE_SHUFFLING;

            update_progression(worker);

        END_JOBCOUNTER_EXCLUSIVE

        /* Generate shuffled read queue from error profile */
        read_queue_set = generate_read_queue_set(&randgen,
                                                 worker->error_profile);
        if (read_queue_set == NULL) {
            LOCK_AND_UPDATE_STATUS(worker, WORKER_STATUS_ERROR);
            fclose(fp);
            return NULL;
        }

        /* Run simulated sequencing and get profiles */
        LOCK_AND_UPDATE_STATUS(worker, WORKER_STATUS_SIMULATING);

        gzrewind(fp);
        simulate_sequencing(fp, read_queue_set, trees);

        free_read_queue_set(read_queue_set);

        /* Increase the job counter for jobs done */
        BEGIN_JOBCOUNTER_EXCLUSIVE

            jobcounter->done++;
            worker->status = WORKER_STATUS_NOT_RUNNING;

            gettimeofday(&tv, NULL);
            worker->last_consumed = TIMEVAL_DIFF(worker->last_started, tv);
            update_progression(worker);

        END_JOBCOUNTER_EXCLUSIVE
    }

    gzclose(fp);

    return (void *)trees;
}

static int
write_permutation_result(const char *prefix, const char *method,
                         CLIPSTATS_ARRAY *arr)
{
    char filename[PATH_MAX];
    struct clipstats *arrptr, *arrend;
    gzFile outf;

    snprintf(filename, PATH_MAX, "%s%s.perm.gz", prefix, method);
    outf = gzopen(filename, "w7");
    if (outf == NULL)
        return -1;

    arrptr = arr->nodes;
    arrend = arr->nodes + arr->end;
    for (; arrptr < arrend; arrptr++)
        gzprintf(outf, "%u %c %f %lu\n", arrptr->depth,
                 NUCLEOBASES[arrptr->basetype], arrptr->value, arrptr->count);

    gzclose(outf);

    return 0;
}

static uint64_t
get_total_read_bases(ERROR_PROFILE *prof)
{
    uint64_t totalbases=0;
    uint32_t pos, base1, base2;

    for (pos = 0; pos < prof->readlength; pos++)
        for (base1 = 0; base1 < NUMBASES; base1++)
            for (base2 = 0; base2 < NUMBASES; base2++)
                totalbases += prof->readdist[pos][base1][base2];

    return totalbases;
}

#ifdef __gnu_linux__
static uint64_t
get_total_memory_size(void)
{
    FILE *fp;
    char label[BUFSIZ], memunit[BUFSIZ];
    uint64_t amount;

    fp = fopen("/proc/meminfo", "r");
    if (fp == NULL)
        return 0;

    if (fscanf(fp, "%s %lu %s", label, &amount, memunit) < 3) {
        fclose(fp);
        return 0;
    }
    fclose(fp);

    if (strcmp(label, "MemTotal:") != 0 || strcmp(memunit, "kB") != 0) {
        fprintf(stderr, "Unexpected first line\n");
        return 0;
    }

    return amount * 1024;
}
#endif

static void
adjust_threads_for_memory(uint64_t totalbases, int *nthreads)
{
#ifdef __gnu_linux__
    uint64_t memsize;
    int recommended_threads;

    memsize = get_total_memory_size();
    if (memsize == 0) {
        printf("Failed to get the total size of memory in /proc/meminfo\n");
        return;
    }

    recommended_threads = (int)(memsize / 1.7 / totalbases);
    if (recommended_threads < 1) {
        printf("Memory is desperately insufficient for this data set. "
               "The program proceeds to run, however it will extremely slow "
               "due to excessive swapping.\n");
        *nthreads = 1;
    }
    else if (recommended_threads < *nthreads) {
        printf("Due to memory insufficiency, number of threads is adjusted "
               "from %d to %d.\n", *nthreads, recommended_threads);
        *nthreads = recommended_threads;
    }
#else
    return;
#endif
}

int
main(int argc, char *argv[])
{
    const char *readpool_input_path, *profile_input_path, *output_prefix;
    ERROR_PROFILE *error_profile=NULL;
    JOB_COUNTER jobcounter;
    int nthreads, i, c, iterations;
    uint64_t coverage_gb=500;
    uint64_t totalbases;
    int quiet=0;

    readpool_input_path = profile_input_path = output_prefix = NULL;
    nthreads = 8;
    iterations = 1000;

    while ((c = getopt(argc, argv, "i:e:m:o:t:d:qh")) != -1)
        switch (c) {
        case 'i':
            readpool_input_path = optarg;
            break;
        case 'e':
            profile_input_path = optarg;
            break;
        case 't':
            nthreads = atoi(optarg);
            break;
        case 'o':
            output_prefix = optarg;
            break;
        case 'd':
            coverage_gb = atol(optarg);
            break;
        case 'q':
            quiet = 1;
            break;
        default:
            usage(argv[0]);
            return 1;
        }

    if (readpool_input_path == NULL) {
        fprintf(stderr, "Input readpool file must be given.\n\n");
        usage(argv[0]);
        return 1;
    }

    if (profile_input_path == NULL) {
        fprintf(stderr, "Input error profile array must be given.\n\n");
        usage(argv[0]);
        return 1;
    }

    if (nthreads < 0 || nthreads > 255) {
        fprintf(stderr, "Number of threads invalid (%d).\n", nthreads);
        return 1;
    }

    if (coverage_gb < 1 || coverage_gb > 100000) {
        fprintf(stderr, "Minimum depth invalid (%lu).\n", coverage_gb);
        return 1;
    }

    if (output_prefix == NULL) {
        fprintf(stderr, "Output files prefix must be given.\n\n");
        usage(argv[0]);
        return 1;
    }

    error_profile = load_error_profile(profile_input_path);
    if (error_profile == NULL)
        goto onError;

    totalbases = get_total_read_bases(error_profile);
    adjust_threads_for_memory(totalbases, &nthreads);
    iterations = (int)ceil(coverage_gb * 1000000000. / totalbases);
    printf("Depth of single run: %lu\n", totalbases);
    printf("Will permutate runs for %d iterations to generate "
        "%.2f Gb.\n\n", iterations, totalbases / 1000000000. * iterations);

    if (iterations < 1 || iterations > 10000) {
        fprintf(stderr, "Number of iterations is determined abnormally: %d\n", iterations);
        return 1;
    }

    jobcounter.total = iterations;
    jobcounter.done = jobcounter.queued = 0;
    jobcounter.nthreads = nthreads;
    gettimeofday(&jobcounter.started, NULL);

    /* main work loop */
    {
        WORKER *workers;
        CLIPSTATS_ARRAY *result_del, *result_mod;
        CLIPSTATS_ARRAY *result_moddel, *result_entropy, *result_t2crate;

        workers = malloc(sizeof(WORKER) * nthreads);
        if (workers == NULL)
            goto onError;

        /* Allocate buffers for collecting the result */
        result_del = clipstatsarray_new(32768);
        result_mod = clipstatsarray_new(32768);
        result_moddel = clipstatsarray_new(32768);
        result_entropy = clipstatsarray_new(32768);
        result_t2crate = clipstatsarray_new(32768);
        if (result_del == NULL || result_mod == NULL || result_moddel == NULL ||
                result_entropy == NULL || result_t2crate == NULL) {
            if (result_del != NULL) clipstatsarray_destroy(result_del);
            if (result_entropy != NULL) clipstatsarray_destroy(result_entropy);
            if (result_moddel != NULL) clipstatsarray_destroy(result_moddel);
            if (result_mod != NULL) clipstatsarray_destroy(result_mod);
            if (result_t2crate != NULL) clipstatsarray_destroy(result_t2crate);
            goto onError;
        }

        /* Initialize per-thread data */
        for (i = 0; i < nthreads; i++) {
            workers[i].workers = workers;
            workers[i].error_profile = error_profile;
            workers[i].readpool_input_path = readpool_input_path;
            workers[i].jobcounter = &jobcounter;
            workers[i].threadid = i;
            workers[i].last_consumed = -1;
            workers[i].status = WORKER_STATUS_NOT_RUNNING;
            workers[i].show_progression = !quiet;
        }

        pthread_mutex_init(&jobcounter.lock, NULL);

        /* Run the thread workers */
        for (i = 0; i < nthreads; i++) {
            if (pthread_create(&workers[i].thread, NULL, &run_clip_permutation,
                               &workers[i]) != 0) {
                fprintf(stderr, "Failed to create thread %d\n", i+1);
                /* TODO: destroy and free threads */
                goto onError;
            }
        }

        /* workers running */

        for (i = 0; i < nthreads; i++) {
            CLIPSTATS_ARRAY *merged_del, *merged_mod;
            CLIPSTATS_ARRAY *merged_moddel, *merged_entropy, *merged_t2crate;
            CLIPSTATS_TREES *trees;

            if (pthread_join(workers[i].thread, (void **)&trees) != 0) {
                fprintf(stderr, "Failed to joining thread %d\n", i+1);
                /* TODO: destroy and free threads */
                goto onError;
            }

            if (trees == NULL) {
                fprintf(stderr, "Thread %d failed to run.\n", i+1);
                /* TODO: destroy and free threads */
                goto onError;
            }

            LOCK_AND_UPDATE_STATUS(&workers[i], WORKER_STATUS_MERGING_OUTPUT);

            merged_del = clipstatsarray_mergetree(result_del, &trees->del);
            merged_mod = clipstatsarray_mergetree(result_mod, &trees->mod);
            merged_moddel = clipstatsarray_mergetree(result_moddel, &trees->moddel);
            merged_entropy = clipstatsarray_mergetree(result_entropy, &trees->entropy);
            merged_t2crate = clipstatsarray_mergetree(result_t2crate, &trees->t2crate);
            if (merged_del == NULL || merged_mod == NULL || merged_t2crate == NULL ||
                    merged_moddel == NULL || merged_entropy == NULL) {
                if (merged_del != NULL) clipstatsarray_destroy(merged_del);
                if (merged_entropy != NULL) clipstatsarray_destroy(merged_entropy);
                if (merged_moddel != NULL) clipstatsarray_destroy(merged_moddel);
                if (merged_mod != NULL) clipstatsarray_destroy(merged_mod);
                if (merged_t2crate!= NULL) clipstatsarray_destroy(merged_t2crate);
                /* TODO: destroy result buffers too. */
                goto onError;
            }

            clipstatsarray_destroy(result_del);
            clipstatsarray_destroy(result_mod);
            clipstatsarray_destroy(result_moddel);
            clipstatsarray_destroy(result_entropy);
            clipstatsarray_destroy(result_t2crate);

            result_del = merged_del; result_mod = merged_mod;
            result_moddel = merged_moddel; result_entropy = merged_entropy;
            result_t2crate = merged_t2crate;

            clipstatstrees_destroy(trees);

            LOCK_AND_UPDATE_STATUS(&workers[i], WORKER_STATUS_FINISHED);
        }

        printf("\n\nWriting out permutation result files.\n");

        write_permutation_result(output_prefix, "del", result_del);
        write_permutation_result(output_prefix, "mod", result_mod);
        write_permutation_result(output_prefix, "moddel", result_moddel);
        write_permutation_result(output_prefix, "entropy", result_entropy);
        write_permutation_result(output_prefix, "t2crate", result_t2crate);

        clipstatsarray_destroy(result_del);
        clipstatsarray_destroy(result_mod);
        clipstatsarray_destroy(result_moddel);
        clipstatsarray_destroy(result_entropy);
        clipstatsarray_destroy(result_t2crate);

        pthread_mutex_destroy(&jobcounter.lock);

        free(workers);

        {
            double elapsed_time;
            struct timeval tv;

            gettimeofday(&tv, NULL);
            elapsed_time = TIMEVAL_DIFF(jobcounter.started, tv);
            printf("\n\n%d iterations successfully finished. "
                   "(%.3lf megabases/sec)\n",
                   iterations,
                   totalbases / 1000000. * iterations / elapsed_time);
        }
    }

    free(error_profile);

    return 0;

  onError:
    if (error_profile != NULL)
        free(error_profile);

    return -1;
}
