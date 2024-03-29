
/* ------------------------------------------------------------
 * This is the file "basic.c" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

#include "h2basic.h"

#include <stdio.h>
#ifdef WIN32
#include <Windows.h>
#include <MMSystem.h>
#pragma comment(lib, "winmm")
#else
#include <sys/times.h>
#include <unistd.h>
#endif

#ifdef USE_CAIRO
#include <cairo/cairo.h>
#include <cairo/cairo-pdf.h>
#endif

#ifdef USE_OPENMP
#include <omp.h>
#endif

#ifdef USE_FREEGLUT
#include <GL/freeglut.h>
#endif

#ifdef USE_GTK3
#include <gtk/gtk.h>
#endif

int       max_pardepth = 0;

/* ------------------------------------------------------------
 * Set up the library
 * ------------------------------------------------------------ */

void
init_h2lib(int *argc, char ***argv)
{
  (void) argc;
  (void) argv;

#ifdef USE_OPENMP
  char     *env;
  int       i, j;

  if (omp_in_parallel()) {
    fprintf(stderr, "Calling init_h2lib in parallel section is forbidden.\n");
    abort();
  }
  omp_set_nested(1);

  i = omp_get_num_procs();
  j = 0;
  while (i > 1) {
    i /= 2;
    j++;
  }
  max_pardepth = j + 1;

  env = getenv("H2_PARDEPTH");
  if (env)
    sscanf(env, "%d", &max_pardepth);
#else
  max_pardepth = 0;
#endif

#ifdef USE_FREEGLUT
  glutInit(argc, *argv);
#endif

#ifdef USE_GTK3
  gtk_disable_setlocale();
  gtk_init(argc, argv);
#endif
}

void
uninit_h2lib()
{

}

/* ------------------------------------------------------------
 * Memory management
 * ------------------------------------------------------------ */

void     *
_h2_allocmem(size_t sz, const char *filename, int line)
{
  void     *ptr;

  ptr = malloc(sz);
  if (ptr == NULL && sz > 0) {
    (void) fprintf(stderr, "Memory allocation of %lu bytes failed in %s:%d\n",
		   (unsigned long) sz, filename, line);
    abort();
  }

  return ptr;
}

uint     *
_h2_allocuint(size_t sz, const char *filename, int line)
{
  uint     *ptr;
  size_t    dsz;

  dsz = sizeof(uint) * sz;
  if (dsz / sizeof(uint) != sz) {
    (void) fprintf(stderr, "Integer overflow in vector allocation in %s:%d\n",
		   filename, line);
    abort();
  }

  ptr = (uint *) malloc(dsz);
  if (ptr == NULL && dsz > 0) {
    (void) fprintf(stderr,
		   "Vector allocation of %lu entries failed in %s:%d\n",
		   (unsigned long) sz, filename, line);
    abort();
  }

  return ptr;
}

double     *
_h2_allocreal(size_t sz, const char *filename, int line)
{
  double     *ptr;
  size_t    dsz;

  dsz = sizeof(double) * sz;
  if (dsz / sizeof(double) != sz) {
    (void) fprintf(stderr, "Integer overflow in vector allocation in %s:%d\n",
		   filename, line);
    abort();
  }

  ptr = (double *) malloc(dsz);
  if (ptr == NULL && dsz > 0) {
    (void) fprintf(stderr,
		   "Vector allocation of %lu entries failed in %s:%d\n",
		   (unsigned long) sz, filename, line);
    abort();
  }

  return ptr;
}

field    *
_h2_allocfield(size_t sz, const char *filename, int line)
{
  field    *ptr;
  size_t    dsz;

  dsz = sizeof(field) * sz;
  if (dsz / sizeof(field) != sz) {
    (void) fprintf(stderr, "Integer overflow in vector allocation in %s:%d\n",
		   filename, line);
    abort();
  }

  ptr = (field *) malloc(dsz);
  if (ptr == NULL && dsz > 0) {
    (void) fprintf(stderr,
		   "Vector allocation of %lu entries failed in %s:%d\n",
		   (unsigned long) sz, filename, line);
    abort();
  }

  return ptr;
}

field    *
_h2_allocmatrix(size_t rows, size_t cols, const char *filename, int line)
{
  field    *ptr;
  size_t    dsz;

  dsz = sizeof(field) * rows * cols;
  if (dsz / sizeof(field) != rows * cols) {
    (void) fprintf(stderr, "Integer overflow in matrix allocation in %s:%d\n",
		   filename, line);
    abort();
  }

  ptr = (field *) malloc(dsz);
  if (ptr == NULL && dsz > 0) {
    (void) fprintf(stderr,
		   "Matrix allocation with %lu rows and %lu columns failed in %s:%d\n",
		   (unsigned long) rows, (unsigned long) cols, filename,
		   line);
    abort();
  }

  return ptr;
}

void
freemem(void *ptr)
{
  free(ptr);
}

/* ------------------------------------------------------------
 * Sorting
 * ------------------------------------------------------------ */

static void
heap_down(int root, int n, bool leq(int, int, void *),
	  void swap(int, int, void *), void *data)
{
  int      child;

  child = 2 * root + 1;
  while (child < n) {
    /* Find larger child */
    if (child + 1 < n && leq(child, child + 1, data))
      child++;

    /* Check heap property */
    if (leq(root, child, data)) {
      swap(root, child, data);
      root = child;
      child = 2 * root + 1;
    }
    else
      break;
  }
}

void
_h2_heapsort(int n, bool leq(int, int, void *),
	     void swap(int, int, void *), void *data)
{
  int      root;

  if (n < 2)
    return;

  /* Build heap */
  root = n / 2;
  while (root > 0) {
    root--;

    heap_down(root, n, leq, swap, data);
  }

  /* Sort */
  while (n > 1) {
    swap(0, n - 1, data);
    n--;
    heap_down(0, n, leq, swap, data);
  }
}

/* ------------------------------------------------------------
 * Timing
 * ------------------------------------------------------------ */

struct _stopwatch {
#ifdef WIN32
  DWORD     start;
  DWORD     current;
#else
#ifdef USE_OPENMP
  double      start;
  double      current;
#else
  struct tms start;
  struct tms current;
#endif
  double      clk_tck;
#endif
};

pstopwatch
new_stopwatch()
{
  pstopwatch sw;

  sw = (pstopwatch) allocmem(sizeof(stopwatch));
#ifndef WIN32
  sw->clk_tck = sysconf(_SC_CLK_TCK);
#endif

  return sw;
}

void
del_stopwatch(pstopwatch sw)
{
  freemem(sw);
}

void
start_stopwatch(pstopwatch sw)
{
#ifdef WIN32
  sw->start = timeGetTime();
#else
#ifdef USE_OPENMP
  sw->start = omp_get_wtime();
#else
  times(&sw->start);
#endif
#endif
}

double
stop_stopwatch(pstopwatch sw)
{
#ifdef WIN32
  sw->current = timeGetTime();
  return (sw->current - sw->start) * 0.001;
#else
#ifdef USE_OPENMP
  sw->current = omp_get_wtime();
  return (sw->current - sw->start);
#else
  times(&sw->current);
  return (sw->current.tms_utime + sw->current.tms_stime - sw->start.tms_utime
	  - sw->start.tms_stime) / sw->clk_tck;
#endif
#endif
}

/* ------------------------------------------------------------
 * Drawing
 * ------------------------------------------------------------ */

#ifdef USE_CAIRO
cairo_t  *
new_cairopdf(const char *filename, double width, double height)
{
  cairo_t  *cr;
  cairo_surface_t *pdf;

  pdf = cairo_pdf_surface_create(filename, width, height);
  cr = cairo_create(pdf);
  cairo_scale(cr, width, height);
  cairo_set_line_width(cr, 0.5 / (width > height ? width : height));
  cairo_surface_destroy(pdf);

  return cr;
}
#endif
