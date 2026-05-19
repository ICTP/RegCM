/*
 * Small pthread shim for RegCM async netCDF output.
 */

#include <pthread.h>
#include <stdint.h>

extern void regcm_async_netcdf_worker(void);

static pthread_t worker_thread;
static int worker_started = 0;

static pthread_mutex_t queue_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t work_cond = PTHREAD_COND_INITIALIZER;
static pthread_cond_t space_cond = PTHREAD_COND_INITIALIZER;
static pthread_cond_t idle_cond = PTHREAD_COND_INITIALIZER;

/* State-serialized netCDF gate; lets foreground callers declare priority. */
static pthread_mutex_t netcdf_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t netcdf_cond = PTHREAD_COND_INITIALIZER;
static int netcdf_worker_active = 0;
static int netcdf_main_active = 0;
static int netcdf_main_waiting = 0;

static void *regcm_async_thread_entry(void *arg)
{
  (void)arg;
  regcm_async_netcdf_worker();
  return 0;
}

int regcm_async_thread_start(void)
{
  int rc;

  if (worker_started) {
    return 0;
  }
  rc = pthread_create(&worker_thread, 0, regcm_async_thread_entry, 0);
  if (rc == 0) {
    worker_started = 1;
  }
  return rc;
}

int regcm_async_thread_join(void)
{
  int rc;

  if (!worker_started) {
    return 0;
  }
  rc = pthread_join(worker_thread, 0);
  if (rc == 0) {
    worker_started = 0;
  }
  return rc;
}

void regcm_async_queue_lock(void)
{
  (void)pthread_mutex_lock(&queue_mutex);
}

void regcm_async_queue_unlock(void)
{
  (void)pthread_mutex_unlock(&queue_mutex);
}

void regcm_async_wait_work(void)
{
  (void)pthread_cond_wait(&work_cond, &queue_mutex);
}

void regcm_async_signal_work(void)
{
  (void)pthread_cond_signal(&work_cond);
}

void regcm_async_broadcast_work(void)
{
  (void)pthread_cond_broadcast(&work_cond);
}

void regcm_async_wait_space(void)
{
  (void)pthread_cond_wait(&space_cond, &queue_mutex);
}

void regcm_async_broadcast_space(void)
{
  (void)pthread_cond_broadcast(&space_cond);
}

void regcm_async_wait_idle(void)
{
  (void)pthread_cond_wait(&idle_cond, &queue_mutex);
}

void regcm_async_broadcast_idle(void)
{
  (void)pthread_cond_broadcast(&idle_cond);
}

void regcm_async_netcdf_lock(void)
{
  (void)pthread_mutex_lock(&netcdf_mutex);
  netcdf_main_waiting++;
  while (netcdf_worker_active || netcdf_main_active) {
    (void)pthread_cond_wait(&netcdf_cond, &netcdf_mutex);
  }
  netcdf_main_waiting--;
  netcdf_main_active = 1;
  (void)pthread_mutex_unlock(&netcdf_mutex);
}

void regcm_async_netcdf_unlock(void)
{
  (void)pthread_mutex_lock(&netcdf_mutex);
  netcdf_main_active = 0;
  (void)pthread_cond_broadcast(&netcdf_cond);
  (void)pthread_mutex_unlock(&netcdf_mutex);
}

void regcm_async_worker_netcdf_lock(void)
{
  (void)pthread_mutex_lock(&netcdf_mutex);
  while (netcdf_main_waiting > 0 || netcdf_main_active || netcdf_worker_active) {
    (void)pthread_cond_wait(&netcdf_cond, &netcdf_mutex);
  }
  netcdf_worker_active = 1;
  (void)pthread_mutex_unlock(&netcdf_mutex);
}

void regcm_async_worker_netcdf_unlock(void)
{
  (void)pthread_mutex_lock(&netcdf_mutex);
  netcdf_worker_active = 0;
  (void)pthread_cond_broadcast(&netcdf_cond);
  (void)pthread_mutex_unlock(&netcdf_mutex);
}

void *regcm_async_ptr_offset(void *base, int64_t offset)
{
  return (void *)((char *)base + offset);
}
