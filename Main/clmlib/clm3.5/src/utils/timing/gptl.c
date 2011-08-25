#include <stdlib.h>        /* malloc */
#include <sys/time.h>      /* gettimeofday */
#include <sys/times.h>     /* times */
#include <unistd.h>        /* gettimeofday */
#include <stdio.h>
#include <string.h>        /* memset, strcmp (via STRMATCH) */
#include <assert.h>

#ifdef UNICOSMP
#include <intrinsics.h>    /* rtc */
#endif

#ifdef SPMD
#include "mpi.h"
#endif

#include "private.h"

/* Use hash table assist to linked list because it is faster */
#define HASH

static Timer **timers = NULL;       /* linked list of timers */
static Timer **lastopen = NULL;     /* last element in list */
static int *max_depth;           /* maximum indentation level encountered */
static int *max_name_len;        /* max length of timer name */

typedef struct {
  unsigned int depth;            /* depth in calling tree */
  int padding[31];               /* padding is to mitigate false cache sharing */
} Nofalse; 
static Nofalse *current_depth;

static int nthreads    = -1;     /* num threads. Init to bad value */
static int maxthreads  = -1;     /* max threads (=nthreads for OMP). Init to bad value */
static int depthlimit  = 99999;  /* max depth for timers (99999 is effectively infinite) */
static bool disabled = false;    /* Timers disabled? */
static bool initialized = false; /* GPTLinitialize has been called */

typedef struct {
  const Option option;           /* wall, cpu, etc. */
  const char *str;               /* descriptive string for printing */
  bool enabled;                  /* flag */
} Settings;

/* Options, print strings, and default enable flags */

static Settings cpustats =      {GPTLcpu,      "Usr       sys       usr+sys   ", false};
static Settings wallstats =     {GPTLwall,     "   Wallclock    max       min          ", true };
static Settings overheadstats = {GPTLoverhead, "Overhead  ",                     true };

static const int tablesize = 128*MAX_CHARS;  /* 128 is size of ASCII char set */
static Hashentry **hashtable;    /* table of entries hashed by sum of chars */
#ifdef DIAG
static unsigned int *novfl;      /* microsecond overflow counter (only when DIAG set) */
#endif
static long ticks_per_sec = -1;  /* clock ticks per second */

#ifdef UNICOSMP
static long long ticks_per_secI = -1;
#endif

/* Local function prototypes */

static void printstats (const Timer *, FILE *, const int, const bool);
static void add (Timer *, const Timer *);
static int get_cpustamp (long *, long *);
static Timer *getentry (const Hashentry *, const char *, int *);
static void printstats_recurse (const Timer *,FILE *,const int,const bool);    
float overhead_sum(const Timer *timer);

/*
** GPTLsetoption: set option value to true or false.
**
** Input arguments:
**   option: option to be set
**   val:    value to which option should be set (for boolean values, 
**           nonzero=true, zero=false)
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTLsetoption (const int option,  /* option */
		   const int val)     /* whether to enable */
{

  if (initialized)
    return GPTLerror ("GPTLsetoption: must be called BEFORE GPTLinitialize\n");

  if (option == GPTLabort_on_error) {
    GPTLset_abort_on_error ((bool) val);
#ifdef DEBUG
    printf ("GPTLsetoption: setting abort on error flag to %d\n", val);
#endif
    return 0;
  }

  switch (option) {
  case GPTLcpu:      
    cpustats.enabled = (bool) val; 
#ifdef DEBUG
    printf ("GPTLsetoption: set cpustats to %d\n", val);
#endif
    return 0;
  case GPTLwall:     
    wallstats.enabled = (bool) val; 
#ifdef DEBUG
    printf ("GPTLsetoption: set wallstats to %d\n", val);
#endif
    return 0;
  case GPTLoverhead: 
    overheadstats.enabled = (bool) val; 
#ifdef DEBUG
    printf ("GPTLsetoption: set overheadstats to %d\n", val);
#endif
    return 0;
  case GPTLdepthlimit: 
    depthlimit = val; 
#ifdef DEBUG
    printf ("GPTLsetoption: set depthlimit to %d\n", val);
#endif
    return 0;
  default:
    break;
  }

#ifdef HAVE_PAPI
  if (GPTL_PAPIsetoption (option, val) == 0)
    return 0;
#endif
  return GPTLerror ("GPTLsetoption: option %d not found\n", option);
}

/*
** GPTLinitialize (): Initialization routine must be called from single-threaded
**   region before any other timing routines may be called.  The need for this
**   routine could be eliminated if not targetting timing library for threaded
**   capability. 
**
** return value: 0 (success) or GPTLerror (failure)
*/

int GPTLinitialize (void)
{
  int i;          /* loop index */
  int t;          /* thread index */
#ifdef UNICOSMP
  extern long long rtc_rate_();
#endif

  if (initialized)
    return GPTLerror ("GPTLinitialize: has already been called\n");

  if (threadinit (&nthreads, &maxthreads) < 0)
    return GPTLerror ("GPTLinitialize: bad return from threadinit\n");

  if (get_thread_num (&nthreads, &maxthreads) > 0) 
    return GPTLerror ("GPTLinitialize: must only be called by master thread\n");

#ifdef UNICOSMP
  ticks_per_secI = rtc_rate_();
#endif

  /* Allocate space for global arrays */

  timers        = (Timer **)     GPTLallocate (maxthreads * sizeof (Timer *));
  lastopen      = (Timer **)     GPTLallocate (maxthreads * sizeof (Timer *));
  current_depth = (Nofalse *)    GPTLallocate (maxthreads * sizeof (Nofalse));
  max_depth     = (int *)        GPTLallocate (maxthreads * sizeof (int));
  max_name_len  = (int *)        GPTLallocate (maxthreads * sizeof (int));
  hashtable     = (Hashentry **) GPTLallocate (maxthreads * sizeof (Hashentry *));
#ifdef DIAG
  novfl         = (unsigned int *) GPTLallocate (maxthreads * sizeof (unsigned int));
#endif

  /* Initialize array values */
#ifdef DEBUG
  printf("GPTL: maxthreads = %d\n", maxthreads);
#endif

  for (t = 0; t < maxthreads; t++) {
    timers[t] = NULL;
    current_depth[t].depth = NULL;
    max_depth[t]     = NULL;
    max_name_len[t]  = NULL;
    hashtable[t] = (Hashentry *) GPTLallocate (tablesize * sizeof (Hashentry));
#ifdef DIAG
    novfl[t] = NULL;
#endif
    for (i = 0; i < tablesize; i++) {
      hashtable[t][i].nument = 0;
      hashtable[t][i].entries = 0;
    }
  }

#ifdef HAVE_PAPI
  if (GPTL_PAPIinitialize (maxthreads) < 0)
    return GPTLerror ("GPTLinitialize: GPTL_PAPIinitialize failure\n");
#endif

  initialized = true;
  return 0;
}

void freelist(Timer *ptr)
{
  if (ptr){
    if(ptr->child)
      freelist(ptr->child);
    if(ptr->next)
      freelist(ptr->next);
    free(ptr);
  }
}

/*
** GPTLfinalize (): Finalization routine must be called from single-threaded
**   region. Free all malloc'd space
**
** return value: 0 (success) or GPTLerror (failure)
*/

int GPTLfinalize (void)
{
  int t;                /* thread index */
  Timer *ptr, *ptrnext; /* ll indices */

  if ( ! initialized)
    return GPTLerror ("GPTLfinalize: GPTLinitialize() has not been called\n");

  if (get_thread_num (&nthreads, &maxthreads) > 0) 
    return GPTLerror ("GPTLfinalize: must only be called by master thread\n");

  for (t = 0; t < maxthreads; ++t) {
    free (hashtable[t]);
    freelist(timers[t]);
  }

  free (timers);
  free (current_depth);
  free (max_depth);
  free (max_name_len);
  free (hashtable);
#ifdef DIAG
  free (novfl);
#endif

  threadfinalize ();
#ifdef HAVE_PAPI
  GPTL_PAPIfinalize (maxthreads);
#endif
  initialized = false;
  return 0;
}

Timer *lastatlevel(Timer *ptr){
#ifdef DEBUG
  printf("GPTL %s \n",ptr->name);
#endif
  if(ptr->next)
    return(lastatlevel(ptr->next));
  else
    return(ptr);
}

Timer *findbyname(Timer *ptr,const char *name, const int t, int *hash_indx){
#ifdef HASH
  ptr = getentry (hashtable[t], name, hash_indx);
  assert (*hash_indx < tablesize);
  return(ptr);
#else
  if(ptr==NULL)
    return(ptr);
  if(STRMATCH(name,ptr->name))
    return(ptr);
  if(ptr->child)
    return(findbyname(ptr->child,name));
  if(ptr->next)
    return(findbyname(ptr->next,name));
  return(NULL);
#endif
}

/*
** GPTLstart: start a timer
**
** Input arguments:
**   name: timer name
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTLstart (const char *name)       /* timer name */
{
#ifdef SPMD
  double wtime1, wtime2;        /* returned from MPI_Wtime() */
#else
  struct timeval tp1, tp2;      /* argument returned from gettimeofday */
#endif

  Timer *ptr=NULL, *ptr2=NULL;            /* linked list pointer */
  Timer **eptr;                 /* for realloc */
  int hash_indx;


  int nchars;                   /* number of characters in timer */
  int t;                        /* thread index (of this thread) */
  int nument;                   /* number of entries for a hash collision */

#ifdef UNICOSMP
  long long nticks1, nticks2;        /* returned from rtc() */

#ifndef SSP
  if (__streaming() == 0) return 0;  /* timers don't work in this situation so disable */
#endif
#endif

  if (disabled){
    return 0;
  }  

  if ((t = get_thread_num (&nthreads, &maxthreads)) < 0)
    return GPTLerror ("GPTLstart\n");

  /*
  ** If current depth exceeds a user-specified limit for print, just
  ** increment and return
  */

  if (current_depth[t].depth >= depthlimit) {
    ++current_depth[t].depth;
    return 0;
  }

  /* 1st calls to overheadstart and gettimeofday are solely for overhead timing */

  if (overheadstats.enabled) {
#ifdef HAVE_PAPI
    (void) GPTL_PAPIoverheadstart (t);
#endif

    if (wallstats.enabled)
#ifdef UNICOSMP
      nticks1 = _rtc();
#else
#ifdef SPMD
      wtime1 = MPI_Wtime();
#else
      gettimeofday (&tp1, 0);
#endif
#endif
  }

  if ( ! initialized)
    return GPTLerror ("GPTLstart: GPTLinitialize has not been called\n");

  /* Look for the requested timer in the current list. */
  
  ptr = findbyname(timers[t],name, t, &hash_indx);

  if (ptr && ptr->onflg)
    return GPTLerror ("GPTLstart thread %d: timer %s was already on: "
		      "not restarting.\n", t, ptr->name);

  ++current_depth[t].depth;
  if (current_depth[t].depth > max_depth[t])
    max_depth[t] = current_depth[t].depth;

  if (ptr==NULL) {

    /* Add a new entry and initialize */

    ptr = (Timer *) GPTLallocate (sizeof (Timer));
    memset (ptr, 0, sizeof (Timer));

    /* Truncate input name if longer than MAX_CHARS characters  */

    nchars = MIN (strlen (name), MAX_CHARS);
    max_name_len[t] = MAX (nchars, max_name_len[t]);

    strncpy (ptr->name, name, nchars);
    ptr->name[nchars] = '\0';
    ptr->depth = current_depth[t].depth;
#ifdef DEBUG
    printf ("GPTL: t=%d timers[t]=%x %s lastopen=%x %x\n",t,timers[t],name,lastopen[t],ptr);
#endif
    if (timers[t]){
      if(lastopen[t] && lastopen[t]->child)
	ptr2 = lastatlevel(lastopen[t]->child);
      else
	ptr2 = NULL;

      if(ptr2){
	ptr->prev=ptr2;
	ptr->next=ptr2->next;
	ptr2->next=ptr;
      }else if(lastopen[t]){
	lastopen[t]->child=ptr;
        ptr->prev=lastopen[t];
      }
    }else{
      timers[t] = ptr;
    }
#ifdef HASH
    ++hashtable[t][hash_indx].nument;
    nument = hashtable[t][hash_indx].nument;

    eptr = realloc (hashtable[t][hash_indx].entries, nument * sizeof (Timer *));
    if ( ! eptr)
      return GPTLerror ("GPTLstart: realloc error\n");

    hashtable[t][hash_indx].entries = eptr;						 
    hashtable[t][hash_indx].entries[nument-1] = ptr;
#endif
  }

  ptr->onflg = true;
  lastopen[t] = ptr;

  if (get_cpustamp (&ptr->cpu.last_utime, &ptr->cpu.last_stime) < 0)
    return GPTLerror ("GPTLstart: get_cpustamp error");
  
  /*
  ** The 2nd system timer call is used both for overhead estimation and
  ** the input timer.  Under UNICOSMP, "ticks" are recorded by rtc instead 
  ** of sec as returned by MPI_Wtime or sec and usec as returned by gettimeofday(). 
  */
  
  if (wallstats.enabled) {
#ifdef UNICOSMP
    nticks2               = _rtc ();
    ptr->wall.last_nticks = nticks2; /* just store rtc to be used and converted later */
    if (overheadstats.enabled)
      ptr->wall.overhead += (float) (nticks2 - nticks1) / (float) ticks_per_secI;
#else
#ifdef SPMD
    wtime2               = MPI_Wtime();
    ptr->wall.last_sec   = wtime2; 
    if (overheadstats.enabled)
      ptr->wall.overhead += wtime2 - wtime1;
#else
    gettimeofday (&tp2, 0);
    ptr->wall.last_sec  = tp2.tv_sec;
    ptr->wall.last_usec = tp2.tv_usec;
    if (overheadstats.enabled)
      ptr->wall.overhead +=       (tp2.tv_sec  - tp1.tv_sec) + 
	                    1.e-6*(tp2.tv_usec - tp1.tv_usec);
#endif
#endif
  }

#ifdef HAVE_PAPI
  if (GPTL_PAPIstart (t, &ptr->aux) < 0)
    return GPTLerror ("GPTLstart: error from GPTL_PAPIstart\n");

  if (overheadstats.enabled)
    (void) GPTL_PAPIoverheadstop (t, &ptr->aux);
#endif

  return (0);
}


/*
** GPTLstop: stop a timer
**
** Input arguments:
**   name: timer name
**
** Return value: 0 (success) or -1 (failure)
*/

int GPTLstop (const char *name) /* timer name */
{
#ifdef SPMD
  double wtime1, wtime2;    /* returned from MPI_Wtime() */
#else
  struct timeval tp1, tp2;  /* argument to gettimeofday() */
  long delta_wtime_sec;     /* wallclock sec change fm GPTLstart() to GPTLstop() */    
  long delta_wtime_usec;    /* wallclock usec change fm GPTLstart() to GPTLstop() */
#endif
  float delta_wtime;        /* floating point wallclock change */
   
  Timer *ptr;               /* linked list pointer */

  int t;                    /* thread number for this process */
  int indx;                 /* index into hash table */

  long usr;                 /* user time (returned from get_cpustamp) */
  long sys;                 /* system time (returned from get_cpustamp) */

#ifdef UNICOSMP
  long long nticks1, nticks2;        /* returned from rtc() */
  long long delta_nticks;            /* change in tick count */

#ifndef SSP
  if (__streaming() == 0) return 0;  /* timers don't work in this situation so disable */
#endif
#endif

  if (disabled){
    return 0;
  }  

  if ((t = get_thread_num (&nthreads, &maxthreads)) < 0)
    return GPTLerror ("GPTLstop\n");

  /*
  ** If current depth exceeds a user-specified limit for print, just
  ** decrement and return
  */

  if (current_depth[t].depth > depthlimit) {
    --current_depth[t].depth;
    return 0;
  }

#ifdef HAVE_PAPI
  if (overheadstats.enabled)
    (void) GPTL_PAPIoverheadstart (t);
#endif

  /*
  ** The 1st system timer call is used both for overhead estimation and
  ** the input timer
  */
    
  if (wallstats.enabled)
#ifdef UNICOSMP
    nticks1 = _rtc ();
#else
#ifdef SPMD
    wtime1 = MPI_Wtime();
#else
    gettimeofday (&tp1, 0);
#endif
#endif

  if (get_cpustamp (&usr, &sys) < 0)
    return GPTLerror (0);

  if ( ! initialized)
    return GPTLerror ("GPTLstop: GPTLinitialize has not been called\n");
  ptr = lastopen[t];


  if(! STRMATCH (name, ptr->name)) /* not closing last opened timer */
    ptr = findbyname(timers[t], name, t, &indx);

  /* reset lastopen pointer to next level */
  for(lastopen[t]=ptr;lastopen[t] && !lastopen[t]->onflg;lastopen[t]=lastopen[t]->prev);

  if ( ! ptr) 
    return GPTLerror ("GPTLstop: thread %d timer for %s had not been started.\n", t, name);

  if ( ! ptr->onflg )
    return GPTLerror ("GPTLstop: thread %d timer %s was already off.\n",t,ptr->name);

#ifdef HAVE_PAPI
  if (GPTL_PAPIstop (t, &ptr->aux) < 0)
    return GPTLerror ("GPTLstart: error from GPTL_PAPIstop\n");
#endif

  --current_depth[t].depth;

  ptr->onflg = false;
  ptr->count++;

  if (wallstats.enabled) {

    /*
    ** Under UNICOSMP, rtc() is being used instead of gettimeofday().
    ** Under SPMD, MPI_Wtime is being used instead of gettimeofday().
    ** In both cases the usec component is not needed.
    */

#ifdef UNICOSMP
    delta_nticks            = nticks1 - ptr->wall.last_nticks;
    delta_wtime             = (float) delta_nticks / (float) ticks_per_secI;
    ptr->wall.accum_nticks += delta_nticks;
#else
#ifdef SPMD
    delta_wtime             = wtime1 - ptr->wall.last_sec;
    ptr->wall.accum_sec    += delta_wtime;
    ptr->wall.accum_usec    = 0;
#else
    delta_wtime_sec         = tp1.tv_sec  - ptr->wall.last_sec;
    delta_wtime_usec        = tp1.tv_usec - ptr->wall.last_usec;
    delta_wtime             = delta_wtime_sec + 1.e-6*delta_wtime_usec;
    ptr->wall.accum_sec    += delta_wtime_sec;
    ptr->wall.accum_usec   += delta_wtime_usec;
#endif
#endif

    if (ptr->count == 1) {
      ptr->wall.max = delta_wtime;
      ptr->wall.min = delta_wtime;
    } else {
      ptr->wall.max = MAX (ptr->wall.max, delta_wtime);
      ptr->wall.min = MIN (ptr->wall.min, delta_wtime);
    }

    /*
    ** Adjust accumulated wallclock values to guard against overflow in the
    ** microsecond accumulator.
    */

    if (ptr->wall.accum_usec > 10000000) {
      ptr->wall.accum_sec  += 10;
      ptr->wall.accum_usec -= 10000000;
#ifdef DIAG
      ++novfl[t];
#endif
    } else if (ptr->wall.accum_usec < -10000000) {
      ptr->wall.accum_sec  -= 10;
      ptr->wall.accum_usec += 10000000;
#ifdef DIAG
      ++novfl[t];
#endif
    }

    /* 2nd system timer call is solely for overhead timing */

    if (overheadstats.enabled) {
#ifdef UNICOSMP
      nticks2 = _rtc ();
      ptr->wall.overhead += (float) (nticks2 - nticks1) / (float) ticks_per_secI;
#else
#ifdef SPMD
      wtime2 = MPI_Wtime();
      ptr->wall.overhead += (wtime2 - wtime1);
#else
      gettimeofday (&tp2, 0);
      ptr->wall.overhead +=       (tp2.tv_sec  - tp1.tv_sec) + 
	                    1.e-6*(tp2.tv_usec - tp1.tv_usec);
#endif
#endif
    }
  }

  if (cpustats.enabled) {
    ptr->cpu.accum_utime += usr - ptr->cpu.last_utime;
    ptr->cpu.accum_stime += sys - ptr->cpu.last_stime;
    ptr->cpu.last_utime   = usr;
    ptr->cpu.last_stime   = sys;
  }

#ifdef HAVE_PAPI
  if (overheadstats.enabled)
    (void) GPTL_PAPIoverheadstop (t, &ptr->aux);
#endif

  return 0;
}

/*
** GPTLenable: enable timers
**
** Return value: 0 (success) 
*/

int GPTLenable (void)
{
  disabled = false;
  return (0);
}

/*
** GPTLdisable: disable timers
**
** Return value: 0 (success) 
*/

int GPTLdisable (void)
{
  disabled = true;
  return (0);
}

/*
** GPTLstamp: Compute timestamp of usr, sys, and wallclock time (seconds)
**
** Output arguments:
**   wall: wallclock
**   usr:  user time
**   sys:  system time
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTLstamp (double *wall, double *usr, double *sys)
{
#ifndef SPMD
  struct timeval tp;         /* argument to gettimeofday */
#endif
  struct tms buf;            /* argument to times */

  if ( ! initialized)
    return GPTLerror ("GPTLstamp: GPTLinitialize has not been called\n");

  *usr = 0;
  *sys = 0;

#ifndef CATAMOUNT
  if (cpustats.enabled){

    if (times (&buf) == -1)
      return GPTLerror ("GPTLstamp: times() failed. Results bogus\n");

    if (ticks_per_sec == -1){
      if ((ticks_per_sec = sysconf (_SC_CLK_TCK)) == -1)
        return GPTLerror ("GPTLinitialize: sysconf (_SC_CLK_TCK) failed\n");
    }

    *usr = buf.tms_utime / (double) ticks_per_sec;
    *sys = buf.tms_stime / (double) ticks_per_sec;
  }
#endif

#ifdef SPMD
  *wall = MPI_Wtime();
#else
  gettimeofday (&tp, 0);
  *wall = tp.tv_sec + 1.e-6*tp.tv_usec;
#endif

  return 0;
}


void reset_all(Timer *ptr)
{
  if (ptr){
    ptr->onflg = false;
    ptr->count = 0;
    memset (&ptr->wall, 0, sizeof (ptr->wall));
    memset (&ptr->cpu, 0, sizeof (ptr->cpu));
#ifdef HAVE_PAPI
    memset (&ptr->aux, 0, sizeof (ptr->aux));
#endif
    if(ptr->child)
      reset_all(ptr->child);
    if(ptr->next)
      reset_all(ptr->next);
  }
}

/*
** GPTLreset: reset all known timers to 0
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTLreset (void)
{
  int t;             /* index over threads */
  Timer *ptr;        /* linked list index */

  if ( ! initialized)
    return GPTLerror ("GPTLreset: GPTLinitialize has not been called\n");

  /* Only allow the master thread to reset timers */

  if (get_thread_num (&nthreads, &maxthreads) > 0) 
    return GPTLerror ("GPTLreset: must only be called by master thread\n");

  for (t = 0; t < nthreads; t++) {
    reset_all(timers[t]);
  }
  printf ("GPTLreset: accumulators for all timers set to zero\n");
  return 0;
}

void printbythread(Timer *ptr, FILE *fp){
  bool foundany, first, found;
  Timer sumstats, *tptr;
  int t, indx;
  /* 
  ** To print sum stats, first create a new timer then copy thread 0
  ** stats into it. then sum using "add", and finally print.
  */

  foundany = false;
  first = true;
  sumstats = *ptr;
  for (t = 1; t < nthreads; ++t) {
    found = false;
    tptr = findbyname(timers[t], ptr->name, t, &indx);
    if(tptr){
      /* Only print thread 0 when this timer found for other threads */
      if (first) {
	first = false;
	fprintf (fp, "%3.3d ", 0);
	printstats (ptr, fp, 0, false);
      }
      
      found = true;
      foundany = true;
      fprintf (fp, "%3.3d ", t);
      printstats (tptr, fp, 0, false);
      add (&sumstats, tptr);
    }
  }
  

  if (foundany) {
    fprintf (fp, "SUM ");
    printstats (&sumstats, fp, 0, false);
    fprintf (fp, "\n");
  }
  if(ptr->child)
    printbythread(ptr->child, fp);
  if(ptr->next)
    printbythread(ptr->next, fp);

}


/* 
** GPTLpr: Print values of all timers
**
** Input arguments:
**   mode:  file open mode: 0 (append), 1 (new)
**   fname: timer output filename
**
** Return value: 0 (success) or GPTLerror (failure)
*/

int GPTLpr (const int mode, const char *fname)
{
  FILE *fp;                /* file handle to write to */
  int i, ii, n, t;         /* indices */
  int me, npes;            /* proc. id and total number of processes */
  int label;               /* timer output label */
  bool first;
  char outfile[21];        /* name of output file: timing.xxxx */
  float *sum;              /* sum of overhead values (per thread) */
  float osum;              /* sum of overhead over threads */
#ifdef SPMD
  MPI_Comm comm_c;         /* MPI communicator */
  int tag;                 /* send/recv tag */
  int destid, srcid;       /* process ids */
  int err;                 /* MPI error return */
  int signal=1;           /* signaling buffer */
  MPI_Status status;       /* MPI status return */
#endif

  if ( ! initialized)
    return GPTLerror ("GPTLpr: GPTLinitialize() has not been called\n");

  if (get_thread_num (&nthreads, &maxthreads) > 0) 
    return GPTLerror ("GPTLpr: must only be called by master thread\n");

  if (mode == 0){
    if ( ! (fp = fopen (fname, "a"))) fp = stderr;
  }
  else{
    if ( ! (fp = fopen (fname, "w"))) fp = stderr;
  }

  sum = (float *) GPTLallocate (nthreads * sizeof (float));

  for (t = 0; t < nthreads; ++t) {
    if (t > 0)
      fprintf (fp, "\n");
    fprintf (fp, "Stats for thread %d:\n", t);

    for (n = 0; n < max_depth[t]; ++n)    /* max indent level (depth starts at 1) */
      fprintf (fp, "  ");
    for (n = 0; n < max_name_len[t]; ++n) /* longest timer name */
      fprintf (fp, " ");

    fprintf (fp, "Called   ");

    /* Print strings for enabled timer types */

    if (cpustats.enabled)
      fprintf (fp, "%s", cpustats.str);
    if (wallstats.enabled) {
      fprintf (fp, "%s", wallstats.str);
      if (overheadstats.enabled)
	fprintf (fp, "%s", overheadstats.str);
    }

#ifdef HAVE_PAPI
    GPTL_PAPIprstr (fp);
#endif

    fprintf (fp, "\n");        /* Done with titles, go to next line */

    printstats_recurse(timers[t],fp,t,true);

    /* Sum of overhead across timers is meaningful */

    if (wallstats.enabled && overheadstats.enabled) {      
      sum[t] = overhead_sum(timers[t]);
      fprintf (fp, "Overhead sum = %9.3f wallclock seconds\n", sum[t]);
    }

  }

  /* Print per-name stats for all threads */

  if (nthreads > 1) {
    fprintf (fp, "\nSame stats sorted by timer for threaded regions:\n");
    fprintf (fp, "Thd ");

    for (n = 0; n < max_name_len[0]; ++n) /* longest timer name */
      fprintf (fp, " ");

    fprintf (fp, "Called   ");

    if (cpustats.enabled)
      fprintf (fp, "%s", cpustats.str);
    if (wallstats.enabled) {
      fprintf (fp, "%s", wallstats.str);
      if (overheadstats.enabled)
	fprintf (fp, "%s", overheadstats.str);
    }

#ifdef HAVE_PAPI
    GPTL_PAPIprstr (fp);
#endif

    fprintf (fp, "\n");

    printbythread(timers[0], fp);


    /* Repeat overhead print in loop over threads */

    if (wallstats.enabled && overheadstats.enabled) {
      osum = 0.;
      for (t = 0; t < nthreads; ++t) {
	fprintf (fp, "OVERHEAD.%3.3d (wallclock seconds) = %9.3f\n", t, sum[t]);
	osum += sum[t];
      }
      fprintf (fp, "OVERHEAD.SUM (wallclock seconds) = %9.3f\n", osum);
    }
  }

#ifdef DIAG
  fprintf (fp, "\n");
  for (t = 0; t < nthreads; ++t) 
    fprintf (fp, "novfl[%d]=%d\n", t, novfl[t]);
#endif

  /* Print hash table stats */

#ifdef HASH  
  for (t = 0; t < nthreads; t++) {
    first = true;
    for (i = 0; i < tablesize; i++) {
      int nument = hashtable[t][i].nument;
      if (nument > 1) {
	if (first) {
	  first = false;
	  fprintf (fp, "\nthread %d had some hash collisions:\n", t);
	}
	fprintf (fp, "hashtable[%d][%d] had %d entries:", t, i, nument);
	for (ii = 0; ii < nument; ii++)
	  fprintf (fp, " %s", hashtable[t][i].entries[ii]->name);
	fprintf (fp, "\n");
      }
    }
  }
#endif

  fprintf (fp, "\n\n");

  fclose (fp);
  free (sum);

  return 0;
}

/*
** overhead_sum: accumulate all overhead times associated
**               with an event tree
**
** Return value: accumlated overhead time
*/

float overhead_sum(const Timer *timer)
{
  float sum;
  if (timer){
    sum = timer->wall.overhead;
    if(timer->child)
      sum += overhead_sum(timer->child);
    if(timer->next)
      sum += overhead_sum(timer->next);
  }
  else sum = 0;
  return(sum);
}

/* 
** printstats: print a single timer
**
** Input arguments:
**   timer:    timer for which to print stats
**   fp:       file descriptor to write to
**   t:        thread number
**   doindent: whether to indent
*/

static void printstats (const Timer *timer,     /* timer to print */
			FILE *fp,               /* file descriptor to write to */
			const int t,            /* thread number */
			const bool doindent)    /* whether indenting will be done */
{
  int i;               /* index */
  int indent;          /* index for indenting */
  int extraspace;      /* for padding to length of longest name */
  long ticks_per_sec;  /* returned from sysconf */
  float usr;           /* user time */
  float sys;           /* system time */
  float usrsys;        /* usr + sys */
  float elapse;        /* elapsed time */

  /* Indent to depth of this timer */

  if (doindent)
    for (indent = 0; indent < timer->depth; ++indent)  /* depth starts at 1 */
      fprintf (fp, "  ");

  fprintf (fp, "%s", timer->name);

  /* Pad to length of longest name */

  extraspace = max_name_len[t] - strlen (timer->name);
  for (i = 0; i < extraspace; ++i)
    fprintf (fp, " ");

  /* Pad to max indent level */

  if (doindent)
    for (indent = timer->depth; indent < max_depth[t]; ++indent)
      fprintf (fp, "  ");

  fprintf (fp, "%8ld ", timer->count);

  if (cpustats.enabled) {
    if ((ticks_per_sec = sysconf (_SC_CLK_TCK)) == -1)
      (void) GPTLerror ("printstats: token _SC_CLK_TCK is not defined\n");
    usr = timer->cpu.accum_utime / (float) ticks_per_sec;
    sys = timer->cpu.accum_stime / (float) ticks_per_sec;
    usrsys = usr + sys;
    fprintf (fp, "%9.3f %9.3f %9.3f ", usr, sys, usrsys);
  }

  if (wallstats.enabled) {
#ifdef UNICOSMP
    elapse = timer->wall.accum_nticks / (float) ticks_per_secI;
#else
#ifdef SPMD
    elapse = timer->wall.accum_sec;
#else
    elapse = timer->wall.accum_sec + 1.e-6*timer->wall.accum_usec;
#endif
#endif
    fprintf (fp, "%12.6f %12.6f %12.6f ", elapse, timer->wall.max, timer->wall.min);
    if (overheadstats.enabled)
      fprintf (fp, "%12.6f ", timer->wall.overhead);
  }

#ifdef HAVE_PAPI
  GPTL_PAPIpr (fp, &timer->aux);
#endif

  fprintf (fp, "\n");
}

/* 
** printstats_recurse: print a linked list of timers
**
** Input arguments:
**   timer:     top of timer list to print
**   fp:        file descriptor to write to
**   t:         thread number
**   doindent:  whether to indent
*/

static void printstats_recurse (const Timer *timer,     /* top of timer list to print */
			FILE *fp,               /* file descriptor to write to */
			const int t,            /* thread number */
			const bool doindent)    /* whether indenting will be done */
{
  if (timer){
    printstats(timer, fp, t, doindent);
    if(timer->child)
      printstats_recurse(timer->child,fp,t,doindent);
    if(timer->next)
      printstats_recurse(timer->next,fp,t,doindent);
  }
    
}

/* 
** add: add the contents of tin to tout
**
** Input arguments:
**   tin:  input timer
** Input/output arguments:
**   tout: output timer summed into
*/

static void add (Timer *tout,   
		 const Timer *tin)
{
  if (wallstats.enabled) {
    tout->count           += tin->count;
    tout->wall.accum_sec  += tin->wall.accum_sec;
    tout->wall.accum_usec += tin->wall.accum_usec;
    if (tout->wall.accum_usec > 10000000) {
      tout->wall.accum_sec  += 10;
      tout->wall.accum_usec -= 10000000;
    } else if (tout->wall.accum_usec < -10000000) {
      tout->wall.accum_sec  -= 10;
      tout->wall.accum_usec += 10000000;
    }
      
    tout->wall.max = MAX (tout->wall.max, tin->wall.max);
    tout->wall.min = MIN (tout->wall.min, tin->wall.min);
    if (overheadstats.enabled)
      tout->wall.overhead += tin->wall.overhead;
  }

  if (cpustats.enabled) {
    tout->cpu.accum_utime += tin->cpu.accum_utime;
    tout->cpu.accum_stime += tin->cpu.accum_stime;
  }
#ifdef HAVE_PAPI
  GPTL_PAPIadd (&tout->aux, &tin->aux);
#endif
}

/*
** get_cpustamp: Invoke the proper system timer and return stats.
**
** Output arguments:
**   usr: user time
**   sys: system time
**
** Return value: 0 (success)
*/

static int get_cpustamp (long *usr, long *sys)
{
  struct tms buf;

  *usr = 0;
  *sys = 0;

  if (cpustats.enabled){

    (void) times (&buf);
    *usr = buf.tms_utime;
    *sys = buf.tms_stime;

  }

  return 0;
}

/*
** getentry: find the entry in the hash table and return a pointer to it.
**
** Input args:
**   hashtable: the hashtable (array)
**   name:      string to be hashed on (specifically, summed)
** Output args:
**   indx:      hashtable index
**
** Return value: pointer to the entry, or NULL if not found
*/

static Timer *getentry (const Hashentry *hashtable, /* hash table */
			const char *name,           /* name to hash */
			int *indx)                  /* hash index */
{
  int i;                 /* loop index */
  const char *c = name;  /* pointer to elements of "name" */

  /* Generate the hash value by summing values of the chars in "name" */

  for (*indx = 0; *c; c++)
    *indx += *c;

  /* 
  ** If nument exceeds 1 there was a hash collision and we must search
  ** linearly through an array for a match
  */

  for (i = 0; i < hashtable[*indx].nument; i++)
    if (STRMATCH (name, hashtable[*indx].entries[i]->name))
      return hashtable[*indx].entries[i];

  return 0;
}
