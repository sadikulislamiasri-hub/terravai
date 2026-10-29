/* Priority-flood depression filling with an imposed epsilon gradient
 * (Barnes, Lehman and Mulla, 2014, Computers & Geosciences 62, 117-127).
 * O(n log n); cells on the grid edge or next to NA cells are outlets. */
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>

typedef struct { double z; int i; } node;

static void heap_push(node *h, int *n, double z, int i) {
  int k = (*n)++;
  while (k > 0) {
    int p = (k - 1) / 2;
    if (h[p].z <= z) break;
    h[k] = h[p]; k = p;
  }
  h[k].z = z; h[k].i = i;
}

static node heap_pop(node *h, int *n) {
  node top = h[0], last = h[--(*n)];
  int k = 0;
  for (;;) {
    int c = 2 * k + 1;
    if (c >= *n) break;
    if (c + 1 < *n && h[c + 1].z < h[c].z) c++;
    if (h[c].z >= last.z) break;
    h[k] = h[c]; k = c;
  }
  if (*n > 0) h[k] = last;
  return top;
}

SEXP C_fill_sinks(SEXP zr, SEXP nrr, SEXP ncr, SEXP epsr) {
  int nr = asInteger(nrr), nc = asInteger(ncr), N = nr * nc;
  double eps = asReal(epsr);
  SEXP out = PROTECT(duplicate(zr));
  double *z = REAL(out);
  char *closed = (char *) R_alloc(N, 1);
  node *h = (node *) R_alloc(N, sizeof(node));
  int hn = 0;
  static const int dr[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
  static const int dc[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
  for (int i = 0; i < N; i++) closed[i] = ISNAN(z[i]) ? 1 : 0;
  for (int c = 0; c < nc; c++) for (int r = 0; r < nr; r++) {
    int i = r + c * nr;
    if (closed[i]) continue;
    int edge = (r == 0 || c == 0 || r == nr - 1 || c == nc - 1);
    for (int k = 0; !edge && k < 8; k++) {
      int j = (r + dr[k]) + (c + dc[k]) * nr;
      if (ISNAN(z[j])) edge = 1;
    }
    if (edge) { closed[i] = 1; heap_push(h, &hn, z[i], i); }
  }
  while (hn > 0) {
    node t = heap_pop(h, &hn);
    int r = t.i % nr, c = t.i / nr;
    for (int k = 0; k < 8; k++) {
      int rr = r + dr[k], cc = c + dc[k];
      if (rr < 0 || cc < 0 || rr >= nr || cc >= nc) continue;
      int j = rr + cc * nr;
      if (closed[j]) continue;
      closed[j] = 1;
      if (z[j] < z[t.i] + eps) z[j] = z[t.i] + eps;
      heap_push(h, &hn, z[j], j);
    }
    if ((hn & 0xFFFF) == 0) R_CheckUserInterrupt();
  }
  UNPROTECT(1);
  return out;
}
