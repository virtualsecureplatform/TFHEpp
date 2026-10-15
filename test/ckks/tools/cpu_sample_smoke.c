/* Five threads each burn300ms of CPU time: expect about750 samples at2ms. */
#include <pthread.h>
#include <stdint.h>
#include <time.h>
static void *burn(void *unused)
{
    (void)unused;
    struct timespec start, now;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
    volatile uint64_t value = 1;
    do {
        for (int i = 0; i < 10000; ++i)
            value = value * UINT64_C(6364136223846793005) + 1;
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &now);
    } while ((now.tv_sec - start.tv_sec) * INT64_C(1000000000) +
                 now.tv_nsec - start.tv_nsec < INT64_C(300000000));
    return 0;
}
int main(void)
{
    pthread_t threads[4];
    for (int i = 0; i < 4; ++i)
        if (pthread_create(&threads[i], 0, burn, 0)) return 1;
    burn(0);
    for (int i = 0; i < 4; ++i)
        if (pthread_join(threads[i], 0)) return 1;
    return 0;
}
