/* qp_json_shim.c -- JSON result logging for Quantum Package (qp3).
 *
 * Owns the yyjson document, the path language, the document envelope,
 * the throttle and the writer thread.  Fortran never sees JSON state:
 * it passes a value and an absolute path, and the two entry points of
 * src/json/json.irp.f choose the file and the run kind.
 *
 * The file is rewritten whole, atomically (temp file + fsync + rename),
 * at most once per QP_JSON_MIN_INTERVAL seconds (default 1.0) and once
 * more at exit.  No signal handler is installed, so validity on disk
 * never depends on the process reaching an orderly exit.
 *
 * The complete public API is the block of prototypes below; it is
 * mirrored one-to-one by the bind(C) interfaces of src/json/qp_json.f90.
 * Not named qp_json.c on purpose: object names come from the file stem
 * and qp_json.f90 already produces qp_json.o.
 */

#define _POSIX_C_SOURCE 200809L

#include <errno.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "yyjson.h"

/* ------------------------------------------------------------------ */
/* Public API                                                          */
/* ------------------------------------------------------------------ */

void qp_json_begin(const char *kind, const char *filename,
                   const char *ezfio_file, int master);
void qp_json_add_int(const char *path, long long v);
void qp_json_add_real(const char *path, double v);
void qp_json_add_bool(const char *path, int v);
void qp_json_add_str(const char *path, const char *v);
void qp_json_add_arr_int(const char *path, const long long *v, int n);
void qp_json_add_arr_real(const char *path, const double *v, int n);

/* ------------------------------------------------------------------ */
/* State                                                               */
/* ------------------------------------------------------------------ */

#define QP_JSON_PATH_MAX 4096
#define QP_JSON_KEY_MAX  256
#define QP_JSON_TICK_NS  100000000LL      /* 100 ms */

#define QP_JSON_WRITE_FLAGS \
  (YYJSON_WRITE_PRETTY_TWO_SPACES | YYJSON_WRITE_NEWLINE_AT_END | \
   YYJSON_WRITE_INF_AND_NAN_AS_NULL)

enum {
  T_INT,
  T_REAL,
  T_BOOL,
  T_STR,
  T_ARR_INT,
  T_ARR_REAL
};

/* Lock order, everywhere: io_mu then mu.  The document is only touched
 * under mu, the file only under io_mu, and mu is never held across a
 * read/write/fsync/rename. */
static pthread_mutex_t mu    = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t io_mu = PTHREAD_MUTEX_INITIALIZER;

static yyjson_mut_doc *doc;                  /* NULL until first add       */
static char     out_path[QP_JSON_PATH_MAX];  /* empty until qp_json_begin  */
static char     run_kind[64];
static char     ezfio_file[QP_JSON_PATH_MAX];
static char     host[256];
static int      enabled = 1;                 /* 0 on non-master ranks      */
static int      dirty;                       /* document newer than file   */
static int      stop_thread;
static int      thread_started;
static int      atexit_done;
static double   min_interval = -1.0;         /* <0: not read from env yet  */
static double   t_last_write = -1e18;        /* CLOCK_MONOTONIC seconds    */
static struct timespec t_start;
static pthread_t writer_tid;

/* ------------------------------------------------------------------ */
/* Small helpers                                                       */
/* ------------------------------------------------------------------ */

static void warn_msg(const char *msg, const char *detail)
{
  if (detail != NULL && detail[0] != '\0')
    fprintf(stderr, "qp_json: %s (%s)\n", msg, detail);
  else
    fprintf(stderr, "qp_json: %s\n", msg);
}

static double now_monotonic(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static double since_start(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)(ts.tv_sec - t_start.tv_sec)
       + 1e-9 * (double)(ts.tv_nsec - t_start.tv_nsec);
}

static void utc_now(char *buf, size_t n)
{
  struct timespec ts;
  struct tm tm;
  time_t sec;

  clock_gettime(CLOCK_REALTIME, &ts);
  sec = ts.tv_sec;
  gmtime_r(&sec, &tm);
  strftime(buf, n, "%Y-%m-%dT%H:%M:%SZ", &tm);
}

/* Copies at most n-1 bytes and always terminates. */
static void copy_str(char *dst, size_t n, const char *src)
{
  if (n == 0) return;
  if (src == NULL) {
    dst[0] = '\0';
    return;
  }
  strncpy(dst, src, n - 1);
  dst[n - 1] = '\0';
}

/* ------------------------------------------------------------------ */
/* Document building.  All of this runs with mu held.                  */
/* ------------------------------------------------------------------ */

static yyjson_mut_doc *ensure_doc(void)
{
  if (doc == NULL) {
    doc = yyjson_mut_doc_new(NULL);
    if (doc == NULL) return NULL;
    yyjson_mut_doc_set_root(doc, yyjson_mut_obj(doc));
  }
  return doc;
}

static yyjson_mut_val *make_value(int type, long long iv, double dv,
                                  const char *sv, const void *data, int n)
{
  yyjson_mut_val *arr;
  const long long *ai;
  const double *ad;
  int i;

  switch (type) {
  case T_INT:
    return yyjson_mut_int(doc, (int64_t)iv);
  case T_REAL:
    return yyjson_mut_real(doc, dv);
  case T_BOOL:
    return yyjson_mut_bool(doc, iv != 0);
  case T_STR:
    return yyjson_mut_strcpy(doc, sv == NULL ? "" : sv);
  case T_ARR_INT:
    arr = yyjson_mut_arr(doc);
    if (arr == NULL) return NULL;
    ai = (const long long *)data;
    for (i = 0; i < n; i++)
      (void)yyjson_mut_arr_add_int(doc, arr, (int64_t)ai[i]);
    return arr;
  case T_ARR_REAL:
    arr = yyjson_mut_arr(doc);
    if (arr == NULL) return NULL;
    ad = (const double *)data;
    for (i = 0; i < n; i++)
      (void)yyjson_mut_arr_add_real(doc, arr, ad[i]);
    return arr;
  default:
    return NULL;
  }
}

/* Inserts one value at an absolute path.  mu held, doc created.
 * Returns 0 on success, -1 on a violation (one stderr line printed). */
static int put_locked(const char *path, int type, long long iv, double dv,
                      const char *sv, const void *data, int n)
{
  yyjson_mut_val *obj, *child, *val;
  const char *p, *q;
  char key[QP_JSON_KEY_MAX];
  size_t len;

  obj = yyjson_mut_doc_get_root(doc);
  if (obj == NULL || !yyjson_mut_is_obj(obj)) return -1;

  p = path;
  for (;;) {
    while (*p == '/') p++;
    if (*p == '\0') {
      warn_msg("ignoring a path without a final key", path);
      return -1;
    }
    q = strchr(p, '/');
    len = (q == NULL) ? strlen(p) : (size_t)(q - p);
    if (len == 0 || len >= sizeof key) {
      warn_msg("ignoring a path with an empty or oversized component", path);
      return -1;
    }
    memcpy(key, p, len);
    key[len] = '\0';

    if (q == NULL) {                    /* last component: the key */
      val = make_value(type, iv, dv, sv, data, n);
      if (val == NULL) {
        warn_msg("could not allocate a JSON value", path);
        return -1;
      }
      /* put, not add: re-logging a path replaces its value.  The key
       * string is copied into the document by yyjson_mut_strcpy. */
      if (!yyjson_mut_obj_put(obj, yyjson_mut_strcpy(doc, key), val)) {
        warn_msg("could not insert a JSON value", path);
        return -1;
      }
      return 0;
    }

    child = yyjson_mut_obj_get(obj, key);
    if (child != NULL && !yyjson_mut_is_obj(child)) {
      warn_msg("replacing a non-object node by an object", path);
      child = NULL;
    }
    if (child == NULL) {
      child = yyjson_mut_obj(doc);
      if (child == NULL
          || !yyjson_mut_obj_put(obj, yyjson_mut_strcpy(doc, key), child)) {
        warn_msg("could not allocate an intermediate object", path);
        return -1;
      }
    }
    obj = child;
    p = q + 1;
  }
}

/* Adds an empty object under a root key when the key is absent. */
static void ensure_empty_obj(const char *key)
{
  yyjson_mut_val *root = yyjson_mut_doc_get_root(doc);

  if (root == NULL || !yyjson_mut_is_obj(root)) return;
  if (yyjson_mut_obj_get(root, key) == NULL)
    (void)yyjson_mut_obj_put(root, yyjson_mut_strcpy(doc, key),
                             yyjson_mut_obj(doc));
}

/* run.timestamp and run.elapsed advance once per actual write, never
 * per json_add.  mu held. */
static void stamp_run(void)
{
  char ts[32];

  utc_now(ts, sizeof ts);
  (void)put_locked("/run/timestamp", T_STR, 0, 0.0, ts, NULL, 0);
  (void)put_locked("/run/elapsed", T_REAL, 0, since_start(), NULL, NULL, 0);
}

/* ------------------------------------------------------------------ */
/* File writing.  Runs with io_mu held.                                */
/* ------------------------------------------------------------------ */

/* Returns 0 on success, an errno value otherwise. */
static int write_file(const char *text, size_t len)
{
  char tmp[QP_JSON_PATH_MAX + 64];
  FILE *f;
  int err = 0;

  if (snprintf(tmp, sizeof tmp, "%s.tmp.%ld", out_path, (long)getpid())
      >= (int)sizeof tmp)
    return ENAMETOOLONG;

  f = fopen(tmp, "wb");
  if (f == NULL) return errno;
  if (len > 0 && fwrite(text, 1, len, f) != len) err = errno ? errno : EIO;
  if (err == 0 && fflush(f) != 0) err = errno ? errno : EIO;
  if (err == 0 && fsync(fileno(f)) != 0) err = errno ? errno : EIO;
  if (fclose(f) != 0 && err == 0) err = errno ? errno : EIO;
  if (err != 0) {
    (void)unlink(tmp);
    return err;
  }
  if (rename(tmp, out_path) != 0) {
    err = errno;
    (void)unlink(tmp);
    return err;
  }
  return 0;
}

/* The only writer.  Takes io_mu then mu; never called with mu held. */
static int do_write(int force)
{
  char *text;
  size_t len = 0;
  int err;

  pthread_mutex_lock(&io_mu);
  pthread_mutex_lock(&mu);
  if (enabled == 0 || doc == NULL || out_path[0] == '\0') {
    pthread_mutex_unlock(&mu);
    pthread_mutex_unlock(&io_mu);
    return 0;
  }
  if (!dirty && !force) {
    pthread_mutex_unlock(&mu);
    pthread_mutex_unlock(&io_mu);
    return 0;
  }
  if (!force && now_monotonic() - t_last_write < min_interval) {
    pthread_mutex_unlock(&mu);
    pthread_mutex_unlock(&io_mu);
    return 0;
  }
  /* Clear before the snapshot: a mutation that arrives while the file
   * is being written re-dirties the document and is picked up next
   * tick. */
  dirty = 0;
  stamp_run();
  text = yyjson_mut_write(doc, QP_JSON_WRITE_FLAGS, &len);
  pthread_mutex_unlock(&mu);

  if (text == NULL) {
    pthread_mutex_lock(&mu);
    dirty = 1;
    pthread_mutex_unlock(&mu);
    warn_msg("could not serialise the document", NULL);
    pthread_mutex_unlock(&io_mu);
    return -1;
  }

  err = write_file(text, len);
  free(text);

  pthread_mutex_lock(&mu);
  t_last_write = now_monotonic();
  if (err != 0) dirty = 1;
  pthread_mutex_unlock(&mu);
  pthread_mutex_unlock(&io_mu);

  if (err != 0)
    warn_msg("could not write the JSON file; previous version intact",
             strerror(err));
  return (err == 0) ? 0 : -1;
}

/* ------------------------------------------------------------------ */
/* Mutation entry point                                                */
/* ------------------------------------------------------------------ */

static void put_item(const char *path, int type, long long iv, double dv,
                     const char *sv, const void *data, int n)
{
  size_t plen;

  if (enabled == 0) return;             /* non-master rank: no work */
  if (path == NULL) {
    warn_msg("ignoring a null path", NULL);
    return;
  }
  plen = strlen(path);
  if (plen == 0 || path[0] != '/' || plen >= QP_JSON_PATH_MAX) {
    warn_msg("ignoring a path that is not absolute or is too long", path);
    return;
  }
  if (strchr(path, ' ') != NULL) {
    warn_msg("ignoring a path that contains a space", path);
    return;
  }

  pthread_mutex_lock(&mu);
  if (ensure_doc() == NULL) {
    pthread_mutex_unlock(&mu);
    warn_msg("out of memory", NULL);
    return;
  }
  if (put_locked(path, type, iv, dv, sv, data, n) == 0)
    dirty = 1;
  pthread_mutex_unlock(&mu);

  /* mu is not recursive: do_write is called after the release. */
  if (min_interval == 0.0) (void)do_write(0);
}

/* ------------------------------------------------------------------ */
/* Initialisation, throttle, exit path                                 */
/* ------------------------------------------------------------------ */

/* Accepts only a plain non-negative decimal number; rejects nan, inf,
 * a negative value and trailing garbage without relying on isfinite(),
 * which fast-math folds away. */
static double parse_interval(const char *s)
{
  char *end = NULL;
  double v;

  if (s == NULL) return 1.0;
  while (*s == ' ' || *s == '\t') s++;
  if (!((*s >= '0' && *s <= '9') || *s == '.' || *s == '+')) return 1.0;
  v = strtod(s, &end);
  if (end == s || *end != '\0') return 1.0;
  if (!(v >= 0.0 && v <= 1.0e9)) return 1.0;
  return v;
}

static void *writer_main(void *arg)
{
  (void)arg;
  for (;;) {
    struct timespec tick;
    int stop;

    tick.tv_sec = 0;
    tick.tv_nsec = QP_JSON_TICK_NS;
    nanosleep(&tick, NULL);
    pthread_mutex_lock(&mu);
    stop = stop_thread;
    pthread_mutex_unlock(&mu);
    if (stop) return NULL;
    (void)do_write(0);
  }
}

static void on_exit(void)
{
  pthread_mutex_lock(&mu);
  if (enabled == 0 || doc == NULL) {
    pthread_mutex_unlock(&mu);
    return;
  }
  stop_thread = 1;
  pthread_mutex_unlock(&mu);

  if (thread_started) {
    pthread_join(writer_tid, NULL);     /* at most one 100 ms tick */
    thread_started = 0;
  }
  if (out_path[0] == '\0') {
    fprintf(stderr, "qp_json: qp_json_begin was never called; "
                    "the collected values are discarded\n");
    return;
  }
  (void)do_write(1);                    /* force, ignore the throttle */
}

void qp_json_begin(const char *kind, const char *filename,
                   const char *ezfio_file_in, int master)
{
  char started[32];

  if (master == 0) {
    pthread_mutex_lock(&mu);
    enabled = 0;
    if (doc != NULL) {
      yyjson_mut_doc_free(doc);
      doc = NULL;
    }
    pthread_mutex_unlock(&mu);
    return;
  }

  if (min_interval < 0.0)               /* read the environment once */
    min_interval = parse_interval(getenv("QP_JSON_MIN_INTERVAL"));

  pthread_mutex_lock(&mu);
  if (ensure_doc() == NULL) {
    pthread_mutex_unlock(&mu);
    warn_msg("out of memory", NULL);
    return;
  }

  if (t_start.tv_sec == 0 && t_start.tv_nsec == 0)
    clock_gettime(CLOCK_MONOTONIC, &t_start);
  if (host[0] == '\0') {
    if (gethostname(host, sizeof host - 1) != 0) host[0] = '\0';
    host[sizeof host - 1] = '\0';
  }
  if (ezfio_file[0] == '\0')
    copy_str(ezfio_file, sizeof ezfio_file, ezfio_file_in);
  copy_str(run_kind, sizeof run_kind, kind);
  utc_now(started, sizeof started);

  (void)put_locked("/schema", T_STR, 0, 0.0, "qp3/1", NULL, 0);
  (void)put_locked("/creator/name", T_STR, 0, 0.0, "Quantum Package", NULL, 0);
  (void)put_locked("/creator/version", T_STR, 0, 0.0, "unknown", NULL, 0);
  (void)put_locked("/provenance/ezfio_file", T_STR, 0, 0.0, ezfio_file,
                   NULL, 0);
  (void)put_locked("/provenance/host", T_STR, 0, 0.0, host, NULL, 0);
  (void)put_locked("/provenance/pid", T_INT, (long long)getpid(), 0.0,
                   NULL, NULL, 0);
  (void)put_locked("/provenance/started", T_STR, 0, 0.0, started, NULL, 0);
  (void)put_locked("/units/energy", T_STR, 0, 0.0, "hartree", NULL, 0);
  (void)put_locked("/units/length", T_STR, 0, 0.0, "bohr", NULL, 0);
  (void)put_locked("/units/time", T_STR, 0, 0.0, "second", NULL, 0);
  ensure_empty_obj("input");
  ensure_empty_obj("results");
  (void)put_locked("/run/kind", T_STR, 0, 0.0, run_kind, NULL, 0);

  if (out_path[0] != '\0') {
    if (strcmp(out_path, filename) == 0) {
      pthread_mutex_unlock(&mu);
      return;
    }
    warn_msg("qp_json_begin called twice with different files; "
             "keeping the first one", out_path);
  } else {
    if (filename == NULL || filename[0] == '\0')
      warn_msg("qp_json_begin called without a file name", NULL);
    else
      copy_str(out_path, sizeof out_path, filename);
  }
  pthread_mutex_unlock(&mu);

  if (out_path[0] == '\0') return;

  if (!atexit_done) {
    (void)atexit(on_exit);
    atexit_done = 1;
  }
  if (min_interval > 0.0) {
    if (pthread_create(&writer_tid, NULL, writer_main, NULL) == 0) {
      thread_started = 1;
    } else {
      warn_msg("could not start the writer thread; writing synchronously",
               NULL);
      min_interval = 0.0;
    }
  }
  /* The file exists, and is valid, from the first millisecond. */
  (void)do_write(1);
}

/* ------------------------------------------------------------------ */
/* Public setters                                                      */
/* ------------------------------------------------------------------ */

void qp_json_add_int(const char *path, long long v)
{
  put_item(path, T_INT, v, 0.0, NULL, NULL, 0);
}

void qp_json_add_real(const char *path, double v)
{
  put_item(path, T_REAL, 0, v, NULL, NULL, 0);
}

void qp_json_add_bool(const char *path, int v)
{
  put_item(path, T_BOOL, v, 0.0, NULL, NULL, 0);
}

void qp_json_add_str(const char *path, const char *v)
{
  put_item(path, T_STR, 0, 0.0, v, NULL, 0);
}

void qp_json_add_arr_int(const char *path, const long long *v, int n)
{
  put_item(path, T_ARR_INT, 0, 0.0, NULL, v, n < 0 ? 0 : n);
}

void qp_json_add_arr_real(const char *path, const double *v, int n)
{
  put_item(path, T_ARR_REAL, 0, 0.0, NULL, v, n < 0 ? 0 : n);
}

/* ------------------------------------------------------------------ */
/* End                                                                 */
/* ------------------------------------------------------------------ */
