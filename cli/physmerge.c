/*
 * physmerge -- C command-line implementation of the physmerge R package.
 *
 * Faithful port of R/physical_merge.R (forward scan + collapse + trim),
 * R/read_sumstat.R (column resolution, TEST / chromosome filters, NA drop)
 * and the ID-export half of R/export.R.
 *
 * Default mode is a single streaming pass: memory use is independent of the
 * number of SNPs (only the current line and one pending block are held).
 * Input that is not already position-sorted within a chromosome is rejected
 * unless --sort is given, which loads the records and stable-sorts them the
 * way R's order() would.
 *
 * Build:  cc -O2 -std=c99 -o physmerge physmerge.c -lz
 */
#define _POSIX_C_SOURCE 200809L
#define _CRT_SECURE_NO_WARNINGS      /* MSVC: fopen/strtok are fine as used here */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#ifdef _MSC_VER
#  define strdup _strdup
#endif
#ifndef PHYSMERGE_NO_ZLIB
#include <zlib.h>
#endif

#define PM_VERSION "0.3.0"
#define PM_BUILD   "c-cli-1"

static void die(const char *fmt, ...);

/* Output files this run created.  A fatal error removes them, so a failed run
   never leaves a truncated block table that looks like a complete one. */
static const char *g_out_path = NULL, *g_snp_path = NULL;
static char g_snpdir_file[4096] = {0};

/* ------------------------------------------------------------------ utils */
static void die(const char *fmt, ...) {
    va_list ap; va_start(ap, fmt);
    fputs("physmerge: error: ", stderr);
    vfprintf(stderr, fmt, ap); fputc('\n', stderr);
    va_end(ap);
    if (g_out_path)      remove(g_out_path);
    if (g_snp_path)      remove(g_snp_path);
    if (g_snpdir_file[0]) remove(g_snpdir_file);
    exit(2);
}

/* Refuse to write over the file we are reading: the output is opened for
   truncation while the input is still being streamed. */
static int same_file(const char *a, const char *b) {
    if (!a || !b) return 0;
    if (!strcmp(a, b)) return 1;
#ifndef _MSC_VER
    struct stat sa, sb;
    if (stat(a, &sa) == 0 && stat(b, &sb) == 0)
        return sa.st_dev == sb.st_dev && sa.st_ino == sb.st_ino;
#endif
    return 0;
}
static void *xmalloc(size_t n) { void *p = malloc(n); if (!p) die("out of memory"); return p; }
static void *xrealloc(void *q, size_t n) { void *p = realloc(q, n); if (!p) die("out of memory"); return p; }

typedef struct { char *p; size_t cap; } Sbuf;
static void sset(Sbuf *s, const char *v) {
    if (!v) { if (s->cap) s->p[0] = '\0'; else { s->cap = 8; s->p = xmalloc(8); s->p[0] = '\0'; } return; }
    size_t n = strlen(v) + 1;
    if (n > s->cap) { s->cap = n * 2; s->p = xrealloc(s->p, s->cap); }
    memcpy(s->p, v, n);
}
static const char *sget(const Sbuf *s) { return s->cap ? s->p : ""; }

/* ----------------------------------------------------------------- reader */
typedef struct {
    FILE *fp;
#ifndef PHYSMERGE_NO_ZLIB
    gzFile gz;
#endif
    int is_gz, eof;
    char *buf; size_t cap, len, pos;
} Reader;

static void rd_open(Reader *r, const char *path) {
    memset(r, 0, sizeof(*r));
    r->cap = 1u << 20;
    r->buf = xmalloc(r->cap);
    int gz = 0;
    if (strcmp(path, "-") != 0) {
        FILE *probe = fopen(path, "rb");
        if (!probe) die("cannot open '%s': %s", path, strerror(errno));
        unsigned char m[2] = {0, 0};
        size_t got = fread(m, 1, 2, probe);
        fclose(probe);
        if (got == 2 && m[0] == 0x1f && m[1] == 0x8b) gz = 1;
    }
#ifdef PHYSMERGE_NO_ZLIB
    if (gz) die("this build has no zlib support; gunzip the input first");
#else
    if (gz) {
        r->gz = gzopen(path, "rb");
        if (!r->gz) die("cannot open '%s'", path);
        gzbuffer(r->gz, 1u << 20);
        r->is_gz = 1;
        return;
    }
#endif
    r->fp = (strcmp(path, "-") == 0) ? stdin : fopen(path, "rb");
    if (!r->fp) die("cannot open '%s': %s", path, strerror(errno));
}
static void rd_close(Reader *r) {
#ifndef PHYSMERGE_NO_ZLIB
    if (r->is_gz) { gzclose(r->gz); free(r->buf); return; }
#endif
    if (r->fp && r->fp != stdin) fclose(r->fp);
    free(r->buf);
}
static int rd_fill(Reader *r) {
    if (r->eof) return 0;
    if (r->pos > 0) { memmove(r->buf, r->buf + r->pos, r->len - r->pos); r->len -= r->pos; r->pos = 0; }
    if (r->len + 1 >= r->cap) { r->cap *= 2; r->buf = xrealloc(r->buf, r->cap); }
    size_t room = r->cap - r->len - 1, n;
#ifndef PHYSMERGE_NO_ZLIB
    if (r->is_gz) { int g = gzread(r->gz, r->buf + r->len, (unsigned)room); if (g < 0) die("gzip read failure"); n = (size_t)g; }
    else
#endif
    n = fread(r->buf + r->len, 1, room, r->fp);
    if (n == 0) { r->eof = 1; return 0; }
    r->len += n;
    return 1;
}
/* returns NUL-terminated line without the newline, or NULL at EOF */
static char *rd_line(Reader *r) {
    for (;;) {
        if (r->pos < r->len) {
            char *s = r->buf + r->pos;
            char *nl = memchr(s, '\n', r->len - r->pos);
            if (nl) {
                size_t L = (size_t)(nl - s);
                r->pos += L + 1;
                if (L && s[L - 1] == '\r') L--;
                s[L] = '\0';
                return s;
            }
        }
        if (!rd_fill(r)) {
            if (r->pos < r->len) {
                char *s = r->buf + r->pos;
                size_t L = r->len - r->pos;
                r->pos = r->len;
                if (L && s[L - 1] == '\r') L--;
                s[L] = '\0';
                return s;
            }
            return NULL;
        }
    }
}

/* ------------------------------------------------------------- line split */
/* splits in place; fills fld[0..want-1] (NULL when the line is short) */
static int split_line(char *s, char sep, char **fld, int want) {
    int k = 0;
    for (int i = 0; i < want; i++) fld[i] = NULL;
    char *p = s;
    for (;;) {
        char *e = strchr(p, sep);
        if (k < want) fld[k] = p;
        k++;
        if (!e) break;
        *e = '\0';
        p = e + 1;
    }
    return k;
}

/* --------------------------------------------------------- numeric parser */
static int parse_num(const char *s, double *out) {
    if (!s || !*s) return 0;
    errno = 0;
    char *end;
    double v = strtod(s, &end);
    if (end == s) return 0;
    while (*end == ' ' || *end == '\t') end++;
    if (*end) return 0;
    if (v != v) return 0;               /* NaN, e.g. "nan" */
    *out = v;
    return 1;
}

/* ------------------------------------------------------------------- core */
typedef struct {
    double start, end, rps_bp, rps_val;
    Sbuf id, line;
} Block;

typedef struct {
    /* configuration */
    double sig_th, window;
    int reward_max;      /* 0 = "min" (p-values), 1 = "max" (statistics) */
    int reset_any;       /* 0 = "best", 1 = "any" */
    int have_chrom, have_id, annotate_full;
    FILE *out, *snpf;
    char *snpdir;
    const char *value_name;
    /* running state */
    Sbuf chrom;
    int chrom_set, in_block, has_cur, has_held;
    double steps, sig_this, last_pos;
    Block cur, held;
    long serial;
    long n_blocks_emitted;
} Core;

static int is_sig(const Core *c, double v)          { return c->reward_max ? (v > c->sig_th) : (v < c->sig_th); }
static int is_better(const Core *c, double v, double b) { return c->reward_max ? (v > b)        : (v < b); }

static void fmt_pos(char *dst, size_t n, double v) {
    if (v == floor(v) && fabs(v) < 1e15) snprintf(dst, n, "%.0f", v);
    else snprintf(dst, n, "%.10g", v);
}

static void emit(Core *c, Block *b) {
    char s1[64], s2[64], s3[64];
    fmt_pos(s1, sizeof s1, b->start);
    fmt_pos(s2, sizeof s2, b->end);
    fmt_pos(s3, sizeof s3, b->rps_bp);
    c->serial++;
    c->n_blocks_emitted++;
    fprintf(c->out, "%ld", c->serial);
    if (c->have_chrom) fprintf(c->out, "\t%s", sget(&c->chrom));
    fprintf(c->out, "\t%s\t%s\t%s", s1, s2, s3);
    if (c->have_id) fprintf(c->out, "\t%s", sget(&b->id));
    fprintf(c->out, "\t%.17g", b->rps_val);
    if (c->annotate_full) fprintf(c->out, "\t%s", sget(&b->line));
    fputc('\n', c->out);

    if (c->snpf) fprintf(c->snpf, "%s\n", c->have_id ? sget(&b->id) : s3);
}

static void block_copy(Block *d, const Block *s) {
    d->start = s->start; d->end = s->end; d->rps_bp = s->rps_bp; d->rps_val = s->rps_val;
    sset(&d->id, sget(&s->id));
    sset(&d->line, sget(&s->line));
}

/* collapse pass + trim pass, streamed with a one-block lookahead */
static void stage(Core *c, Block *b) {
    if (c->has_held) {
        if ((b->rps_bp - c->held.rps_bp) < c->window) {           /* .collapse_blocks */
            if (b->end > c->held.end) c->held.end = b->end;
            if (is_better(c, b->rps_val, c->held.rps_val)) {
                c->held.rps_bp = b->rps_bp; c->held.rps_val = b->rps_val;
                sset(&c->held.id, sget(&b->id)); sset(&c->held.line, sget(&b->line));
            }
            return;
        }
        if (c->held.end > b->start) c->held.end = b->start;       /* trim pass */
        emit(c, &c->held);
    }
    block_copy(&c->held, b);
    c->has_held = 1;
}

static void open_snpdir_file(Core *c) {
    if (!c->snpdir) return;
    if (c->snpf) fclose(c->snpf);
    /* the chromosome comes from the input file, so keep it out of the path */
    char safe[64]; size_t k = 0;
    for (const char *q = sget(&c->chrom); *q && k + 1 < sizeof safe; q++)
        safe[k++] = (*q == '/' || *q == '\\' || *q == '.') ? '_' : *q;
    safe[k] = '\0';
    char path[4096];
    snprintf(path, sizeof path, "%s/snp_ch%s.txt", c->snpdir, safe);
    c->snpf = fopen(path, "w");
    if (!c->snpf) die("cannot write '%s': %s", path, strerror(errno));
    snprintf(g_snpdir_file, sizeof g_snpdir_file, "%s", path);
}

static void close_block(Core *c, double last_inblock_pos) {
    c->cur.end = last_inblock_pos + c->steps;
    c->in_block = 0;
    c->steps = c->window;
    c->sig_this = c->sig_th;
    stage(c, &c->cur);
}
static void open_block(Core *c, double pos, double val, const char *id, const char *line) {
    c->cur.start = (pos - c->window) > 0 ? (pos - c->window) : 0;   /* max(0, pos - window) */
    c->cur.end = 0;
    c->cur.rps_bp = pos;
    c->cur.rps_val = val;
    sset(&c->cur.id, id); sset(&c->cur.line, line);
    c->in_block = 1;
    c->steps = c->window;
    c->sig_this = val;
}
static void end_chrom(Core *c) {
    if (!c->chrom_set) return;
    if (c->in_block) close_block(c, c->last_pos);
    if (c->has_held) { emit(c, &c->held); c->has_held = 0; }
    c->in_block = 0;
}

static void core_push(Core *c, const char *chrom, double pos, double val,
                      const char *id, const char *line) {
    if (!c->chrom_set || (c->have_chrom && strcmp(sget(&c->chrom), chrom) != 0)) {
        end_chrom(c);
        sset(&c->chrom, chrom);
        c->chrom_set = 1;
        c->steps = c->window;
        c->sig_this = c->sig_th;
        c->last_pos = pos;
        if (c->snpdir) open_snpdir_file(c);
    }
    if (!c->in_block) {
        if (is_sig(c, val)) open_block(c, pos, val, id, line);
    } else {
        double remaining = c->steps - (pos - c->last_pos);
        if (remaining <= 0) {
            close_block(c, c->last_pos);
            if (is_sig(c, val)) open_block(c, pos, val, id, line);
        } else {
            c->steps = remaining;
            if (c->reset_any && is_sig(c, val)) {
                c->steps = c->window;
                if (is_better(c, val, c->sig_this)) {
                    c->sig_this = val;
                    c->cur.rps_bp = pos; c->cur.rps_val = val;
                    sset(&c->cur.id, id); sset(&c->cur.line, line);
                }
            } else if (is_better(c, val, c->sig_this)) {
                c->sig_this = val;
                c->steps = c->window;
                c->cur.rps_bp = pos; c->cur.rps_val = val;
                sset(&c->cur.id, id); sset(&c->cur.line, line);
            }
        }
    }
    c->last_pos = pos;
}

/* ------------------------------------------------------- --sort buffering */
typedef struct { int chrom_rank; double pos, val; unsigned long idx; size_t id_off, line_off; } Rec;
typedef struct { char *p; size_t len, cap; } Arena;
static size_t arena_put(Arena *a, const char *s) {
    if (!s) return (size_t)-1;
    size_t n = strlen(s) + 1;
    if (a->len + n > a->cap) { while (a->len + n > a->cap) a->cap = a->cap ? a->cap * 2 : (1u << 16); a->p = xrealloc(a->p, a->cap); }
    size_t off = a->len; memcpy(a->p + off, s, n); a->len += n; return off;
}
static const char *arena_get(const Arena *a, size_t off) { return off == (size_t)-1 ? "" : a->p + off; }
static int rec_cmp(const void *A, const void *B) {
    const Rec *a = A, *b = B;
    if (a->chrom_rank != b->chrom_rank) return a->chrom_rank < b->chrom_rank ? -1 : 1;
    if (a->pos != b->pos) return a->pos < b->pos ? -1 : 1;
    return a->idx < b->idx ? -1 : (a->idx > b->idx ? 1 : 0);   /* stable, like order() */
}

/* ------------------------------------------------------------------- main */
static const char *USAGE =
"physmerge " PM_VERSION " (" PM_BUILD ") -- panel-free physical locus merging\n"
"\n"
"Usage: physmerge --input FILE [options]\n"
"\n"
"Input\n"
"  -i, --input FILE       summary statistics (plain, .gz, or '-' for stdin)\n"
"  -f, --format FMT       plink2 (default) | gpcm | custom\n"
"      --chrom-col NAME   chromosome column   (plink2/gpcm default: #CHROM)\n"
"      --pos-col NAME     position column     (default: POS)\n"
"      --id-col NAME      SNP id column       (default: ID; 'NA' to disable)\n"
"      --value-col NAME   value column        (plink2: P, gpcm: P_HPI)\n"
"      --sep CHAR         field separator; default auto-detect from header\n"
"      --no-chrom         ignore the chromosome column entirely\n"
"\n"
"Filters (as in read_sumstat)\n"
"      --test-col NAME    default TEST\n"
"      --test-val VALUE   default ADD\n"
"      --test-filter      force the TEST filter on  (default: on for plink2)\n"
"      --no-test-filter   force the TEST filter off\n"
"      --chrom LIST       comma-separated chromosomes to keep\n"
"\n"
"Merging (as in physical_merge)\n"
"  -s, --sig-th NUM       significance threshold (default 5e-8)\n"
"  -w, --window NUM       window in bp (default 500000)\n"
"  -r, --reward min|max   min for p-values (default), max for statistics\n"
"      --reset-on best|any  window reset rule (default best)\n"
"\n"
"Output\n"
"  -o, --out FILE         block table (default stdout)\n"
"      --snp-list FILE    representative SNP ids, one per line\n"
"      --snp-list-dir DIR one snp_ch<CHR>.txt per chromosome\n"
"      --annotate-full    append the full original line of each lead SNP\n"
"      --no-header        suppress the output header line\n"
"      --sort             buffer and sort the input (needed if unsorted)\n"
"  -q, --quiet            suppress progress messages\n"
"  -h, --help             this help;  --version  print version\n";

int main(int argc, char **argv) {
    const char *in_path = NULL, *out_path = NULL, *snp_path = NULL, *snp_dir = NULL;
    const char *format = "plink2";
    const char *chrom_col = NULL, *pos_col = NULL, *id_col = NULL, *value_col = NULL;
    const char *test_col = "TEST", *test_val = "ADD", *chrom_keep = NULL;
    int test_filter = -1, want_sort = 0, quiet = 0, no_header = 0, annotate_full = 0, no_chrom = 0;
    char sep = 0;
    double sig_th = 5e-8, window = 500000.0;
    int reward_max = 0, reset_any = 0;

#define NEXTARG(name) (++i < argc ? argv[i] : (die("missing value for %s", name), (char*)NULL))
    for (int i = 1; i < argc; i++) {
        char *a = argv[i];
        if (!strcmp(a, "-h") || !strcmp(a, "--help")) { fputs(USAGE, stdout); return 0; }
        else if (!strcmp(a, "--version")) { printf("physmerge %s (%s)\n", PM_VERSION, PM_BUILD); return 0; }
        else if (!strcmp(a, "-i") || !strcmp(a, "--input"))  in_path = NEXTARG("--input");
        else if (!strcmp(a, "-o") || !strcmp(a, "--out"))    out_path = NEXTARG("--out");
        else if (!strcmp(a, "-f") || !strcmp(a, "--format")) format = NEXTARG("--format");
        else if (!strcmp(a, "--chrom-col")) chrom_col = NEXTARG("--chrom-col");
        else if (!strcmp(a, "--pos-col"))   pos_col   = NEXTARG("--pos-col");
        else if (!strcmp(a, "--id-col"))    id_col    = NEXTARG("--id-col");
        else if (!strcmp(a, "--value-col")) value_col = NEXTARG("--value-col");
        else if (!strcmp(a, "--test-col"))  test_col  = NEXTARG("--test-col");
        else if (!strcmp(a, "--test-val"))  test_val  = NEXTARG("--test-val");
        else if (!strcmp(a, "--test-filter"))    test_filter = 1;
        else if (!strcmp(a, "--no-test-filter")) test_filter = 0;
        else if (!strcmp(a, "--chrom"))     chrom_keep = NEXTARG("--chrom");
        else if (!strcmp(a, "--no-chrom"))  no_chrom = 1;
        else if (!strcmp(a, "-s") || !strcmp(a, "--sig-th")) { if (!parse_num(NEXTARG("--sig-th"), &sig_th)) die("--sig-th must be numeric"); }
        else if (!strcmp(a, "-w") || !strcmp(a, "--window")) { if (!parse_num(NEXTARG("--window"), &window)) die("--window must be numeric"); }
        else if (!strcmp(a, "-r") || !strcmp(a, "--reward")) { const char *v = NEXTARG("--reward");
            if (!strcmp(v, "max")) reward_max = 1; else if (!strcmp(v, "min")) reward_max = 0; else die("`reward` must be either 'min' or 'max'"); }
        else if (!strcmp(a, "--reset-on")) { const char *v = NEXTARG("--reset-on");
            if (!strcmp(v, "any")) reset_any = 1; else if (!strcmp(v, "best")) reset_any = 0; else die("`reset_on` must be either 'best' or 'any'"); }
        else if (!strcmp(a, "--snp-list"))     snp_path = NEXTARG("--snp-list");
        else if (!strcmp(a, "--snp-list-dir")) snp_dir  = NEXTARG("--snp-list-dir");
        else if (!strcmp(a, "--annotate-full")) annotate_full = 1;
        else if (!strcmp(a, "--no-header"))    no_header = 1;
        else if (!strcmp(a, "--sort"))         want_sort = 1;
        else if (!strcmp(a, "--sep")) { const char *v = NEXTARG("--sep");
            if (!strcmp(v, "tab") || !strcmp(v, "\\t")) sep = '\t';
            else if (!strcmp(v, "space")) sep = ' ';
            else if (strlen(v) == 1) sep = v[0]; else die("--sep must be a single character"); }
        else if (!strcmp(a, "-q") || !strcmp(a, "--quiet")) quiet = 1;
        else die("unknown option '%s' (try --help)", a);
    }
    if (!in_path) { fputs(USAGE, stderr); return 2; }
    if (window <= 0) die("`window` must be a single positive numeric value.");
    if (snp_path && snp_dir) die("use either --snp-list or --snp-list-dir, not both");

    /* format defaults, mirroring read_sumstat() */
    int dflt_test_filter = 0;
    if (!strcmp(format, "plink2")) {
        if (!chrom_col) chrom_col = "#CHROM";
        if (!pos_col)   pos_col   = "POS";
        if (!id_col)    id_col    = "ID";
        if (!value_col) value_col = "P";
        dflt_test_filter = 1;
    } else if (!strcmp(format, "gpcm")) {
        if (!chrom_col) chrom_col = "#CHROM";
        if (!pos_col)   pos_col   = "POS";
        if (!id_col)    id_col    = "ID";
        if (!value_col) value_col = "P_HPI";
    } else if (!strcmp(format, "custom")) {
        if (!pos_col || !value_col || (!chrom_col && !no_chrom))
            die("For format = 'custom', you must supply --chrom-col, --pos-col and --value-col.");
    } else die("--format must be plink2, gpcm or custom");
    if (test_filter < 0) test_filter = dflt_test_filter;
    if (id_col && !strcmp(id_col, "NA")) id_col = NULL;
    if (no_chrom) chrom_col = NULL;
    /* read_sumstat normalizes #CHROM -> CHROM */
    if (chrom_col && !strcmp(chrom_col, "#CHROM")) chrom_col = "CHROM";

    if (same_file(in_path, out_path))
        die("--out is the same file as --input; choose a different output path");
    if (same_file(in_path, snp_path))
        die("--snp-list is the same file as --input; choose a different output path");

    Reader rd; rd_open(&rd, in_path);
    char *hdr = rd_line(&rd);
    if (!hdr) die("empty input file");
    if (!sep) {
        if (strchr(hdr, '\t')) sep = '\t';
        else if (strchr(hdr, ',')) sep = ',';
        else sep = ' ';
    }
    /* count header fields, then index them */
    int nf = 1; for (char *p = hdr; *p; p++) if (*p == sep) nf++;
    char *hdr_copy = strdup(hdr);          /* split_line() destroys hdr */
    if (!hdr_copy) die("out of memory");
    char **hf = xmalloc((size_t)nf * sizeof(char *));
    split_line(hdr, sep, hf, nf);
    int i_chrom = -1, i_pos = -1, i_id = -1, i_val = -1, i_test = -1;
    for (int k = 0; k < nf; k++) {
        const char *h = hf[k]; if (!h) continue;
        if (!strcmp(h, "#CHROM")) h = "CHROM";                 /* same normalization */
        if (chrom_col && i_chrom < 0 && !strcmp(h, chrom_col)) i_chrom = k;
        if (pos_col   && i_pos   < 0 && !strcmp(h, pos_col))   i_pos   = k;
        if (id_col    && i_id    < 0 && !strcmp(h, id_col))    i_id    = k;
        if (value_col && i_val   < 0 && !strcmp(h, value_col)) i_val   = k;
        if (test_col  && i_test  < 0 && !strcmp(h, test_col))  i_test  = k;
    }
    {
        char miss[512]; miss[0] = '\0';
        if (chrom_col && i_chrom < 0) { strncat(miss, chrom_col, 100); strcat(miss, ", "); }
        if (i_pos < 0) { strncat(miss, pos_col, 100); strcat(miss, ", "); }
        if (i_val < 0) { strncat(miss, value_col, 100); strcat(miss, ", "); }
        if (id_col && i_id < 0) { strncat(miss, id_col, 100); strcat(miss, ", "); }
        if (miss[0]) { size_t L = strlen(miss); miss[L - 2] = '\0'; die("Column(s) not found: %s", miss); }
    }
    if (test_filter && i_test < 0) {
        if (!quiet) fprintf(stderr, "physmerge: warning: test_col '%s' not found; TEST filter skipped.\n", test_col);
        test_filter = 0;
    }
    int max_idx = i_pos; if (i_val > max_idx) max_idx = i_val;
    if (i_chrom > max_idx) max_idx = i_chrom;
    if (i_id > max_idx) max_idx = i_id;
    if (test_filter && i_test > max_idx) max_idx = i_test;
    int want = max_idx + 1;
    char **fld = xmalloc((size_t)want * sizeof(char *));

    /* --chrom keep list */
    char **keep = NULL; int n_keep = 0;
    if (chrom_keep) {
        char *cp = strdup(chrom_keep);
        for (char *tok = strtok(cp, ","); tok; tok = strtok(NULL, ",")) {
            keep = xrealloc(keep, (size_t)(n_keep + 1) * sizeof(char *));
            keep[n_keep++] = strdup(tok);
        }
        free(cp);
    }

    Core c; memset(&c, 0, sizeof c);
    c.sig_th = sig_th; c.window = window; c.reward_max = reward_max; c.reset_any = reset_any;
    c.have_chrom = (i_chrom >= 0); c.have_id = (i_id >= 0); c.annotate_full = annotate_full;
    c.value_name = value_col;
    c.out = out_path ? fopen(out_path, "w") : stdout;
    if (!c.out) die("cannot write '%s': %s", out_path, strerror(errno));
    g_out_path = out_path;
    if (snp_path) { c.snpf = fopen(snp_path, "w"); if (!c.snpf) die("cannot write '%s': %s", snp_path, strerror(errno)); g_snp_path = snp_path; }
    if (snp_dir) c.snpdir = strdup(snp_dir);
    c.steps = window; c.sig_this = sig_th;

    if (!no_header) {
        fputs("serial", c.out);
        if (c.have_chrom) fputs("\tCHROM", c.out);
        fputs("\tstart\tend\trps_BP", c.out);
        if (c.have_id) fputs("\trps_ID", c.out);
        fprintf(c.out, "\trps_%s", value_col);
        if (annotate_full) fprintf(c.out, "\t%s", hdr_copy);   /* original columns */
        fputc('\n', c.out);
    }

    unsigned long n_read = 0, n_kept = 0, n_test_drop = 0, n_chrom_drop = 0, n_na = 0;
    char *linecopy = NULL; size_t linecopy_cap = 0;

    Rec *recs = NULL; size_t n_rec = 0, cap_rec = 0;
    Arena arena; memset(&arena, 0, sizeof arena); arena.p = NULL;
    char **cnames = NULL; int n_cn = 0;

    /* one shared per-record handler */
    Sbuf prev_chrom; memset(&prev_chrom, 0, sizeof prev_chrom);
    int prev_set = 0; double prev_pos = 0;

    char *ln;
    while ((ln = rd_line(&rd)) != NULL) {
        if (!*ln) continue;
        n_read++;
        size_t llen = strlen(ln);
        if (annotate_full) {
            if (llen + 1 > linecopy_cap) { linecopy_cap = (llen + 1) * 2; linecopy = xrealloc(linecopy, linecopy_cap); }
            memcpy(linecopy, ln, llen + 1);
        }
        split_line(ln, sep, fld, want);
        if (test_filter) {
            const char *tv = fld[i_test];
            if (!tv || strcmp(tv, test_val) != 0) { n_test_drop++; continue; }
        }
        const char *ch = (i_chrom >= 0 && fld[i_chrom]) ? fld[i_chrom] : "";
        if (n_keep) {
            int ok = 0;
            for (int k = 0; k < n_keep; k++) if (!strcmp(ch, keep[k])) { ok = 1; break; }
            if (!ok) { n_chrom_drop++; continue; }
        }
        double pos, val;
        if (!parse_num(fld[i_pos], &pos) || !parse_num(fld[i_val], &val)) { n_na++; continue; }
        const char *id = (i_id >= 0 && fld[i_id]) ? fld[i_id] : "";
        n_kept++;

        if (want_sort) {
            int rank = -1;
            for (int k = 0; k < n_cn; k++) if (!strcmp(cnames[k], ch)) { rank = k; break; }
            if (rank < 0) { cnames = xrealloc(cnames, (size_t)(n_cn + 1) * sizeof(char *)); cnames[n_cn] = strdup(ch); rank = n_cn++; }
            if (n_rec == cap_rec) { cap_rec = cap_rec ? cap_rec * 2 : 65536; recs = xrealloc(recs, cap_rec * sizeof(Rec)); }
            recs[n_rec].chrom_rank = rank; recs[n_rec].pos = pos; recs[n_rec].val = val;
            recs[n_rec].idx = (unsigned long)n_rec;
            recs[n_rec].id_off = c.have_id ? arena_put(&arena, id) : (size_t)-1;
            recs[n_rec].line_off = annotate_full ? arena_put(&arena, linecopy) : (size_t)-1;
            n_rec++;
            continue;
        }
        /* streaming: verify the ordering assumption R's order() would enforce */
        if (prev_set) {
            int same = c.have_chrom ? (strcmp(sget(&prev_chrom), ch) == 0) : 1;
            if (same) {
                if (pos < prev_pos)
                    die("input is not position-sorted (chromosome %s: %.0f after %.0f).\n"
                        "       re-run with --sort, or sort the file first.", ch, pos, prev_pos);
            } else {
                for (int k = 0; k < n_cn; k++) if (!strcmp(cnames[k], ch))
                    die("chromosome %s appears in more than one block of the file.\n"
                        "       re-run with --sort, or sort the file first.", ch);
                cnames = xrealloc(cnames, (size_t)(n_cn + 1) * sizeof(char *));
                cnames[n_cn++] = strdup(ch);
            }
        } else {
            cnames = xrealloc(cnames, sizeof(char *)); cnames[0] = strdup(ch); n_cn = 1;
        }
        sset(&prev_chrom, ch); prev_pos = pos; prev_set = 1;

        core_push(&c, ch, pos, val, id, annotate_full ? linecopy : "");
    }
    rd_close(&rd);

    if (want_sort) {
        qsort(recs, n_rec, sizeof(Rec), rec_cmp);
        for (size_t k = 0; k < n_rec; k++)
            core_push(&c, cnames[recs[k].chrom_rank], recs[k].pos, recs[k].val,
                      arena_get(&arena, recs[k].id_off), arena_get(&arena, recs[k].line_off));
    }
    end_chrom(&c);

    if (c.out != stdout) fclose(c.out);
    if (c.snpf) fclose(c.snpf);
    g_out_path = g_snp_path = NULL; g_snpdir_file[0] = '\0';   /* run succeeded */

    if (!quiet) {
        if (test_filter) fprintf(stderr, "physmerge: TEST filter: kept %lu of %lu rows where %s = '%s'.\n",
                                 n_read - n_test_drop, n_read, test_col, test_val);
        if (n_chrom_drop) fprintf(stderr, "physmerge: %lu row(s) dropped by the chromosome filter.\n", n_chrom_drop);
        if (n_na) fprintf(stderr, "physmerge: %lu row(s) dropped (NA in position or value).\n", n_na);
        fprintf(stderr, "physmerge: %lu SNPs -> %ld blocks (window=%.0f, sig_th=%g, reward=%s, reset_on=%s).\n",
                n_kept, c.n_blocks_emitted, window, sig_th, reward_max ? "max" : "min", reset_any ? "any" : "best");
    }
    return 0;
}
