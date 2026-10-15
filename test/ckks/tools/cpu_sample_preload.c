/* Optional Linux/x86-64 CPU sampling diagnostic; no performance counters needed.
 * Build: gcc -std=gnu11 -O2 -shared -fPIC cpu_sample_preload.c -o sampler.so
 *        -ldl -lrt -pthread
 * Run: CKKS_CPU_SAMPLES=/absolute/output-prefix LD_PRELOAD=./sampler.so program
 * Records instruction addresses, counts and module mappings; no operand data.
 * Per-thread CPU timers sample every2ms. Wall-clock buckets are100ms wide.
 * Timed qualification runs should be separate from instrumented profiling runs.
 */
#define _GNU_SOURCE
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <signal.h>
#include <time.h>
#include <ucontext.h>
#include <dlfcn.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <string.h>
#include <errno.h>
#if !defined(__x86_64__)
#error This diagnostic sampler currently records x86-64 RIP.
#endif
_Static_assert(ATOMIC_LLONG_LOCK_FREE == 2, "signal samples require lock-free counters");
#define SLOTS (1u<<21)
struct sample_slot { _Atomic uint64_t key, count; };
static struct sample_slot samples[SLOTS];
static _Atomic uint64_t dropped, created, failed, active;
static _Atomic int stopping;
static _Atomic uint64_t start_ns;
static int sample_signal;
static timer_t main_timer;
static int have_main_timer;
typedef int (*create_function)(pthread_t*,const pthread_attr_t*,void*(*)(void*),void*);
static _Atomic(create_function) next_create;
static create_function resolve_create(void){create_function f=atomic_load(&next_create);if(!f){f=(create_function)dlsym(RTLD_NEXT,"pthread_create");atomic_store(&next_create,f);}return f;}
static uint64_t monotonic_ns(void) {
 struct timespec ts;clock_gettime(CLOCK_MONOTONIC,&ts);
 return (uint64_t)ts.tv_sec*1000000000ull+ts.tv_nsec;
}
static void sample(int sig,siginfo_t*info,void*context) {
 (void)sig;(void)info;
 if(atomic_load_explicit(&stopping,memory_order_relaxed))return;
 const uint64_t pc=(uint64_t)((ucontext_t*)context)->uc_mcontext.gregs[REG_RIP];
 const uint64_t bucket=(monotonic_ns()-atomic_load_explicit(&start_ns,memory_order_relaxed))/100000000ull;
 if(pc>=(1ull<<48)||bucket>=65535){atomic_fetch_add_explicit(&dropped,1,memory_order_relaxed);return;}
 const uint64_t key=(pc<<16)|bucket;
 const uint64_t hash=(key^(key>>33))*0xff51afd7ed558ccdull;
 uint32_t index=(uint32_t)(hash^(hash>>32))&(SLOTS-1);
 for(unsigned probe=0;probe<64;++probe,index=(index+1)&(SLOTS-1)) {
  uint64_t expected=0;
  if(atomic_compare_exchange_strong_explicit(&samples[index].key,&expected,key,memory_order_relaxed,memory_order_relaxed)||expected==key){
   atomic_fetch_add_explicit(&samples[index].count,1,memory_order_relaxed);return;
  }
 }
 atomic_fetch_add_explicit(&dropped,1,memory_order_relaxed);
}
static int install_timer(timer_t*timer) {
 struct sigevent event={0};event.sigev_notify=SIGEV_THREAD_ID;event.sigev_signo=sample_signal;event._sigev_un._tid=(pid_t)syscall(SYS_gettid);
 if(timer_create(CLOCK_THREAD_CPUTIME_ID,&event,timer)){atomic_fetch_add(&failed,1);return 0;}
 struct itimerspec interval={0};interval.it_value.tv_nsec=2000000;interval.it_interval=interval.it_value;
 if(timer_settime(*timer,0,&interval,0)){timer_delete(*timer);atomic_fetch_add(&failed,1);return 0;}
 atomic_fetch_add(&created,1);atomic_fetch_add(&active,1);return 1;
}
struct launch {void*(*function)(void*);void*argument;};
struct cleanup {timer_t timer;int installed;};
static void cleanup_timer(void*argument){struct cleanup*c=argument;if(c->installed){timer_delete(c->timer);atomic_fetch_sub(&active,1);}}
static void*thread_start(void*argument) {
 struct launch launch=*(struct launch*)argument;free(argument);
 struct cleanup cleanup={0};cleanup.installed=install_timer(&cleanup.timer);void*result;
 pthread_cleanup_push(cleanup_timer,&cleanup);
 result=launch.function(launch.argument);
 pthread_cleanup_pop(1);
 return result;
}
int pthread_create(pthread_t*thread,const pthread_attr_t*attributes,void*(*function)(void*),void*argument) {
 create_function create=resolve_create();if(!create)return ENOSYS;
 if(atomic_load(&stopping))return create(thread,attributes,function,argument);
 struct launch*launch=malloc(sizeof(*launch));if(!launch)return ENOMEM;
 launch->function=function;launch->argument=argument;
 int status=create(thread,attributes,thread_start,launch);if(status)free(launch);return status;
}
__attribute__((constructor)) static void begin_sampling(void) {
 if(!getenv("CKKS_CPU_SAMPLES")){atomic_store(&stopping,1);return;}
 (void)resolve_create();sample_signal=SIGRTMIN+6;start_ns=monotonic_ns();
 struct sigaction action={0};action.sa_sigaction=sample;action.sa_flags=SA_SIGINFO|SA_RESTART;sigemptyset(&action.sa_mask);
 if(sigaction(sample_signal,&action,0)){atomic_store(&stopping,1);return;}
 have_main_timer=install_timer(&main_timer);
}
__attribute__((destructor)) static void end_sampling(void) {
 const char*prefix=getenv("CKKS_CPU_SAMPLES");if(!prefix)return;
 atomic_store(&stopping,1);if(have_main_timer){timer_delete(main_timer);atomic_fetch_sub(&active,1);}
 char path[4096];snprintf(path,sizeof(path),"%s.csv",prefix);FILE*out=fopen(path,"w");if(!out)return;
 fprintf(out,"# start_ns=%llu bucket_ns=100000000 interval_cpu_ns=2000000 threads=%llu failures=%llu active=%llu dropped=%llu\n",(unsigned long long)start_ns,(unsigned long long)atomic_load(&created),(unsigned long long)atomic_load(&failed),(unsigned long long)atomic_load(&active),(unsigned long long)atomic_load(&dropped));
 fprintf(out,"bucket,pc,count\n");
 for(size_t i=0;i<SLOTS;++i){uint64_t key=atomic_load(&samples[i].key);if(key)fprintf(out,"%llu,0x%llx,%llu\n",(unsigned long long)(key&65535),(unsigned long long)(key>>16),(unsigned long long)atomic_load(&samples[i].count));}
 fclose(out);snprintf(path,sizeof(path),"%s.maps",prefix);out=fopen(path,"w");FILE*in=fopen("/proc/self/maps","r");
 if(out&&in){char buffer[8192];size_t count;while((count=fread(buffer,1,sizeof(buffer),in)))fwrite(buffer,1,count,out);}
 if(in)fclose(in);
 if(out)fclose(out);
}
