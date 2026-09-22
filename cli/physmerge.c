#define _POSIX_C_SOURCE 200809L
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#ifdef _MSC_VER
# define strdup _strdup
#endif
#ifndef PHYSMERGE_NO_ZLIB
#include <zlib.h>
#endif

#define PM_VERSION "0.3.0"
#define PM_BUILD "c-cli-1"

static void die(const char *fmt, ...);

static const char *g_out = NULL, *g_snp = NULL;
static char g_snpf[4096] = {0};

static void die(const char *fmt, ...) {
    va_list ap; va_start(ap, fmt);
    fputs("physmerge: error: ", stderr);
    vfprintf(stderr, fmt, ap); fputc('\n', stderr);
    va_end(ap);
    if (g_out) remove(g_out);
    if (g_snp) remove(g_snp);
    if (g_snpf[0]) remove(g_snpf);
    exit(2);
}

static int samef(const char *a, const char *b) {
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

static char *fsep(char *p, char sep, int ws) {
    if (!ws) return strchr(p, sep);
    for (; *p; p++) if (*p == ' ' || *p == '\t') return p;
    return NULL;
}

static int split(char *s, char sep, char **fld, int want) {
    int k = 0, ws = (sep == ' ');
    for (int i = 0; i < want; i++) fld[i] = NULL;
    char *p = s;
    if (ws) while (*p == ' ' || *p == '\t') p++;
    for (;;) {
        char *val = p, *rest = p;
        if (*p == '"') {
            char *r = p + 1, *w = p;
            val = w;
            while (*r) {
                if (*r == '"') {
                    if (r[1] == '"') { *w++ = '"'; r += 2; continue; }
                    r++;
                    break;
                }
                *w++ = *r++;
            }
            rest = r;
            *w = '\0';
        }
        char *e = fsep(rest, sep, ws);
        if (k < want) fld[k] = val;
        k++;
        if (!e) break;
        *e = '\0';
        p = e + 1;
        if (ws) while (*p == ' ' || *p == '\t') p++;
    }
    return k;
}

static int pnum(const char *s, double *out) {
    if (!s || !*s) return 0;
    errno = 0;
    char *end;
    double v = strtod(s, &end);
    if (end == s) return 0;
    while (*end == ' ' || *end == '\t') end++;
    if (*end) return 0;
    if (v != v) return 0;
    *out = v;
    return 1;
}

typedef struct {
    double start, end, bp, val;
    Sbuf id, line;
} Block;

typedef struct {

    double sig, window;
    int rmax;
    int rany;
    int havech, haveid, annot;
    FILE *out, *snpf;
    char *snpdir;
    const char *vname;
    char *hdr;

    Sbuf chrom;
    int chset, inblk, hascur, hashld;
    double steps, best, lpos;
    Block cur, held;
    long serial;
    long nblk;
} Core;

static int is_sig(const Core *c, double v) { return c->rmax ? (v > c->sig) : (v < c->sig); }
static int is_better(const Core *c, double v, double b) { return c->rmax ? (v > b) : (v < b); }

static void fmt_pos(char *dst, size_t n, double v) {
    if (v == floor(v) && fabs(v) < 1e15) snprintf(dst, n, "%.0f", v);
    else snprintf(dst, n, "%.10g", v);
}

static void puthdr(Core *c) {
    if (!c->hdr) return;
    fputs(c->hdr, c->out);
    free(c->hdr);
    c->hdr = NULL;
}

static void emit(Core *c, Block *b) {
    char s1[64], s2[64], s3[64];
    puthdr(c);
    if (b->end < b->start) b->end = b->start;
    fmt_pos(s1, sizeof s1, b->start);
    fmt_pos(s2, sizeof s2, b->end);
    fmt_pos(s3, sizeof s3, b->bp);
    c->serial++;
    c->nblk++;
    fprintf(c->out, "%ld", c->serial);
    if (c->havech) fprintf(c->out, "\t%s", sget(&c->chrom));
    fprintf(c->out, "\t%s\t%s\t%s", s1, s2, s3);
    if (c->haveid) fprintf(c->out, "\t%s", sget(&b->id));
    fprintf(c->out, "\t%.10g", b->val);
    if (c->annot) fprintf(c->out, "\t%s", sget(&b->line));
    fputc('\n', c->out);

    if (c->snpf) fprintf(c->snpf, "%s\n", c->haveid ? sget(&b->id) : s3);
}

static void bcopy(Block *d, const Block *s) {
    d->start = s->start; d->end = s->end; d->bp = s->bp; d->val = s->val;
    sset(&d->id, sget(&s->id));
    sset(&d->line, sget(&s->line));
}

static void stage(Core *c, Block *b) {
    if (c->hashld) {
        if ((b->bp - c->held.bp) < c->window) {
            if (b->end > c->held.end) c->held.end = b->end;
            if (is_better(c, b->val, c->held.val)) {
                c->held.bp = b->bp; c->held.val = b->val;
                sset(&c->held.id, sget(&b->id)); sset(&c->held.line, sget(&b->line));
            }
            return;
        }
        if (c->held.end > b->start) c->held.end = b->start;
        emit(c, &c->held);
    }
    bcopy(&c->held, b);
    c->hashld = 1;
}

static void open_snpf(Core *c) {
    if (!c->snpdir) return;
    if (c->snpf) fclose(c->snpf);

    const char *src = sget(&c->chrom);
    char safe[256]; size_t k = 0;
    if (strlen(src) + 1 > sizeof safe)
        die("chromosome name is too long for a file name (%zu characters)", strlen(src));
    for (const char *q = src; *q; q++)
        safe[k++] = (*q == '/' || *q == '\\' || *q == '.') ? '_' : *q;
    safe[k] = '\0';
    char path[4096];
    snprintf(path, sizeof path, "%s/snp_ch%s.txt", c->snpdir, safe);
    c->snpf = fopen(path, "w");
    if (!c->snpf) die("cannot write '%s': %s", path, strerror(errno));
    snprintf(g_snpf, sizeof g_snpf, "%s", path);
}

static void closeb(Core *c, double lastp) {
    c->cur.end = lastp + c->steps;
    c->inblk = 0;
    c->steps = c->window;
    c->best = c->sig;
    stage(c, &c->cur);
}
static void openb(Core *c, double pos, double val, const char *id, const char *line) {
    c->cur.start = (pos - c->window) > 0 ? (pos - c->window) : 0;
    c->cur.end = 0;
    c->cur.bp = pos;
    c->cur.val = val;
    sset(&c->cur.id, id); sset(&c->cur.line, line);
    c->inblk = 1;
    c->steps = c->window;
    c->best = val;
}
static void endch(Core *c) {
    if (!c->chset) return;
    if (c->inblk) closeb(c, c->lpos);
    if (c->hashld) { emit(c, &c->held); c->hashld = 0; }
    c->inblk = 0;
}

static void push(Core *c, const char *chrom, double pos, double val,
                      const char *id, const char *line) {
    if (!c->chset || (c->havech && strcmp(sget(&c->chrom), chrom) != 0)) {
        endch(c);
        sset(&c->chrom, chrom);
        c->chset = 1;
        c->steps = c->window;
        c->best = c->sig;
        c->lpos = pos;
        if (c->snpdir) open_snpf(c);
    }
    if (!c->inblk) {
        if (is_sig(c, val)) openb(c, pos, val, id, line);
    } else {
        double remaining = c->steps - (pos - c->lpos);
        if (remaining <= 0) {
            closeb(c, c->lpos);
            if (is_sig(c, val)) openb(c, pos, val, id, line);
        } else {
            c->steps = remaining;
            if (c->rany && is_sig(c, val)) {
                c->steps = c->window;
                if (is_better(c, val, c->best)) {
                    c->best = val;
                    c->cur.bp = pos; c->cur.val = val;
                    sset(&c->cur.id, id); sset(&c->cur.line, line);
                }
            } else if (is_better(c, val, c->best)) {
                c->best = val;
                c->steps = c->window;
                c->cur.bp = pos; c->cur.val = val;
                sset(&c->cur.id, id); sset(&c->cur.line, line);
            }
        }
    }
    c->lpos = pos;
}

typedef struct { int chrom_rank; double pos, val; unsigned long idx; size_t id_off, line_off; } Rec;
typedef struct { char *p; size_t len, cap; } Arena;
static size_t aput(Arena *a, const char *s) {
    if (!s) return (size_t)-1;
    size_t n = strlen(s) + 1;
    if (a->len + n > a->cap) { while (a->len + n > a->cap) a->cap = a->cap ? a->cap * 2 : (1u << 16); a->p = xrealloc(a->p, a->cap); }
    size_t off = a->len; memcpy(a->p + off, s, n); a->len += n; return off;
}
static const char *aget(const Arena *a, size_t off) { return off == (size_t)-1 ? "" : a->p + off; }
static int rcmp(const void *A, const void *B) {
    const Rec *a = A, *b = B;
    if (a->chrom_rank != b->chrom_rank) return a->chrom_rank < b->chrom_rank ? -1 : 1;
    if (a->pos != b->pos) return a->pos < b->pos ? -1 : 1;
    return a->idx < b->idx ? -1 : (a->idx > b->idx ? 1 : 0);
}

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
"      --reset-on best|any  window reset rule (default any)\n"
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
    const char *inp = NULL, *outp = NULL, *snpp = NULL, *snpd = NULL;
    const char *format = "plink2";
    const char *ccol = NULL, *pcol = NULL, *icol = NULL, *vcol = NULL;
    const char *tcol = "TEST", *tval = "ADD", *chkeep = NULL;
    int tfilt = -1, dosort = 0, quiet = 0, nohdr = 0, annot = 0, nochr = 0;
    char sep = 0;
    double sig = 5e-8, window = 500000.0;
    int rmax = 0, rany = 1;

#define NEXTARG(name) (++i < argc ? argv[i] : (die("missing value for %s", name), (char*)NULL))
    for (int i = 1; i < argc; i++) {
        char *a = argv[i];
        if (!strcmp(a, "-h") || !strcmp(a, "--help")) { fputs(USAGE, stdout); return 0; }
        else if (!strcmp(a, "--version")) { printf("physmerge %s (%s)\n", PM_VERSION, PM_BUILD); return 0; }
        else if (!strcmp(a, "-i") || !strcmp(a, "--input")) inp = NEXTARG("--input");
        else if (!strcmp(a, "-o") || !strcmp(a, "--out")) outp = NEXTARG("--out");
        else if (!strcmp(a, "-f") || !strcmp(a, "--format")) format = NEXTARG("--format");
        else if (!strcmp(a, "--chrom-col")) ccol = NEXTARG("--chrom-col");
        else if (!strcmp(a, "--pos-col")) pcol = NEXTARG("--pos-col");
        else if (!strcmp(a, "--id-col")) icol = NEXTARG("--id-col");
        else if (!strcmp(a, "--value-col")) vcol = NEXTARG("--value-col");
        else if (!strcmp(a, "--test-col")) tcol = NEXTARG("--test-col");
        else if (!strcmp(a, "--test-val")) tval = NEXTARG("--test-val");
        else if (!strcmp(a, "--test-filter")) tfilt = 1;
        else if (!strcmp(a, "--no-test-filter")) tfilt = 0;
        else if (!strcmp(a, "--chrom")) chkeep = NEXTARG("--chrom");
        else if (!strcmp(a, "--no-chrom")) nochr = 1;
        else if (!strcmp(a, "-s") || !strcmp(a, "--sig-th")) { if (!pnum(NEXTARG("--sig-th"), &sig)) die("--sig-th must be numeric"); }
        else if (!strcmp(a, "-w") || !strcmp(a, "--window")) { if (!pnum(NEXTARG("--window"), &window)) die("--window must be numeric"); }
        else if (!strcmp(a, "-r") || !strcmp(a, "--reward")) { const char *v = NEXTARG("--reward");
            if (!strcmp(v, "max")) rmax = 1; else if (!strcmp(v, "min")) rmax = 0; else die("`reward` must be either 'min' or 'max'"); }
        else if (!strcmp(a, "--reset-on")) { const char *v = NEXTARG("--reset-on");
            if (!strcmp(v, "any")) rany = 1; else if (!strcmp(v, "best")) rany = 0; else die("`reset_on` must be either 'best' or 'any'"); }
        else if (!strcmp(a, "--snp-list")) snpp = NEXTARG("--snp-list");
        else if (!strcmp(a, "--snp-list-dir")) snpd = NEXTARG("--snp-list-dir");
        else if (!strcmp(a, "--annotate-full")) annot = 1;
        else if (!strcmp(a, "--no-header")) nohdr = 1;
        else if (!strcmp(a, "--sort")) dosort = 1;
        else if (!strcmp(a, "--sep")) { const char *v = NEXTARG("--sep");
            if (!strcmp(v, "tab") || !strcmp(v, "\\t")) sep = '\t';
            else if (!strcmp(v, "space")) sep = ' ';
            else if (strlen(v) == 1) sep = v[0]; else die("--sep must be a single character"); }
        else if (!strcmp(a, "-q") || !strcmp(a, "--quiet")) quiet = 1;
        else die("unknown option '%s' (try --help)", a);
    }
    if (!inp) { fputs(USAGE, stderr); return 2; }
    if (window <= 0) die("`window` must be a single positive numeric value.");
    if (snpp && snpd) die("use either --snp-list or --snp-list-dir, not both");

    int dtfilt = 0;
    if (!strcmp(format, "plink2")) {
        if (!ccol) ccol = "#CHROM";
        if (!pcol) pcol = "POS";
        if (!icol) icol = "ID";
        if (!vcol) vcol = "P";
        dtfilt = 1;
    } else if (!strcmp(format, "gpcm")) {
        if (!ccol) ccol = "#CHROM";
        if (!pcol) pcol = "POS";
        if (!icol) icol = "ID";
        if (!vcol) vcol = "P_HPI";
    } else if (!strcmp(format, "custom")) {
        if (!pcol || !vcol || (!ccol && !nochr))
            die("For format = 'custom', you must supply --chrom-col, --pos-col and --value-col.");
    } else die("--format must be plink2, gpcm or custom");
    if (tfilt < 0) tfilt = dtfilt;
    if (icol && !strcmp(icol, "NA")) icol = NULL;
    if (nochr) ccol = NULL;

    if (ccol && !strcmp(ccol, "#CHROM")) ccol = "CHROM";

    if (samef(inp, outp))
        die("--out is the same file as --input; choose a different output path");
    if (samef(inp, snpp))
        die("--snp-list is the same file as --input; choose a different output path");

    Reader rd; rd_open(&rd, inp);
    char *hdr = rd_line(&rd);
    if (!hdr) die("empty input file");

    if ((unsigned char)hdr[0] == 0xEF && (unsigned char)hdr[1] == 0xBB &&
        (unsigned char)hdr[2] == 0xBF) hdr += 3;
    if (!sep) {
        if (strchr(hdr, '\t')) sep = '\t';
        else if (strchr(hdr, ',')) sep = ',';
        else sep = ' ';
    }

    int nf = 1; for (char *p = hdr; *p; p++) if (*p == sep) nf++;
    char *hcopy = strdup(hdr);
    if (!hcopy) die("out of memory");
    char **hf = xmalloc((size_t)nf * sizeof(char *));
    split(hdr, sep, hf, nf);
    int i_chrom = -1, i_pos = -1, i_id = -1, i_val = -1, i_test = -1;
    for (int k = 0; k < nf; k++) {
        const char *h = hf[k]; if (!h) continue;
        if (!strcmp(h, "#CHROM")) h = "CHROM";
        if (ccol && i_chrom < 0 && !strcmp(h, ccol)) i_chrom = k;
        if (pcol && i_pos < 0 && !strcmp(h, pcol)) i_pos = k;
        if (icol && i_id < 0 && !strcmp(h, icol)) i_id = k;
        if (vcol && i_val < 0 && !strcmp(h, vcol)) i_val = k;
        if (tcol && i_test < 0 && !strcmp(h, tcol)) i_test = k;
    }
    {
        char miss[512]; miss[0] = '\0';
        if (ccol && i_chrom < 0) { strncat(miss, ccol, 100); strcat(miss, ", "); }
        if (i_pos < 0) { strncat(miss, pcol, 100); strcat(miss, ", "); }
        if (i_val < 0) { strncat(miss, vcol, 100); strcat(miss, ", "); }
        if (icol && i_id < 0) { strncat(miss, icol, 100); strcat(miss, ", "); }
        if (miss[0]) { size_t L = strlen(miss); miss[L - 2] = '\0'; die("Column(s) not found: %s", miss); }
    }
    if (tfilt && i_test < 0) {
        if (!quiet) fprintf(stderr, "physmerge: warning: test_col '%s' not found; TEST filter skipped.\n", tcol);
        tfilt = 0;
    }
    int maxi = i_pos; if (i_val > maxi) maxi = i_val;
    if (i_chrom > maxi) maxi = i_chrom;
    if (i_id > maxi) maxi = i_id;
    if (tfilt && i_test > maxi) maxi = i_test;
    int want = maxi + 1;
    char **fld = xmalloc((size_t)want * sizeof(char *));

    char **keep = NULL; int n_keep = 0;
    if (chkeep) {
        char *cp = strdup(chkeep);
        for (char *tok = strtok(cp, ","); tok; tok = strtok(NULL, ",")) {
            keep = xrealloc(keep, (size_t)(n_keep + 1) * sizeof(char *));
            keep[n_keep++] = strdup(tok);
        }
        free(cp);
    }

    Core c; memset(&c, 0, sizeof c);
    c.sig = sig; c.window = window; c.rmax = rmax; c.rany = rany;
    c.havech = (i_chrom >= 0); c.haveid = (i_id >= 0); c.annot = annot;
    c.vname = vcol;
    c.out = outp ? fopen(outp, "w") : stdout;
    if (!c.out) die("cannot write '%s': %s", outp, strerror(errno));
    g_out = outp;
    if (snpp) { c.snpf = fopen(snpp, "w"); if (!c.snpf) die("cannot write '%s': %s", snpp, strerror(errno)); g_snp = snpp; }
    if (snpd) c.snpdir = strdup(snpd);
    c.steps = window; c.best = sig;

    if (!nohdr) {
        size_t hn = strlen(vcol) + (annot && hcopy ? strlen(hcopy) : 0) + 128;
        char *h = xmalloc(hn);
        h[0] = '\0';
        strcat(h, "serial");
        if (c.havech) strcat(h, "\tCHROM");
        strcat(h, "\tstart\tend\trps_BP");
        if (c.haveid) strcat(h, "\trps_ID");
        strcat(h, "\trps_"); strcat(h, vcol);
        if (annot && hcopy) { strcat(h, "\t"); strcat(h, hcopy); }
        strcat(h, "\n");
        c.hdr = h;
    }

    unsigned long nread = 0, nkept = 0, ntdrop = 0, ncdrop = 0, nna = 0;
    int wneg = 0;
    char *lcopy = NULL; size_t lcap = 0;

    Rec *recs = NULL; size_t n_rec = 0, cap_rec = 0;
    Arena arena; memset(&arena, 0, sizeof arena); arena.p = NULL;
    char **cnames = NULL; int n_cn = 0;

    Sbuf pch; memset(&pch, 0, sizeof pch);
    int pset = 0; double ppos = 0;

    char *ln;
    while ((ln = rd_line(&rd)) != NULL) {
        if (!*ln) continue;
        nread++;
        size_t llen = strlen(ln);
        if (annot) {
            if (llen + 1 > lcap) { lcap = (llen + 1) * 2; lcopy = xrealloc(lcopy, lcap); }
            memcpy(lcopy, ln, llen + 1);
        }
        split(ln, sep, fld, want);
        if (tfilt) {
            const char *tv = fld[i_test];
            if (!tv || strcmp(tv, tval) != 0) { ntdrop++; continue; }
        }
        const char *ch = (i_chrom >= 0 && fld[i_chrom]) ? fld[i_chrom] : "";
        if (n_keep) {
            int ok = 0;
            for (int k = 0; k < n_keep; k++) if (!strcmp(ch, keep[k])) { ok = 1; break; }
            if (!ok) { ncdrop++; continue; }
        }
        double pos, val;
        if (!pnum(fld[i_pos], &pos) || !pnum(fld[i_val], &val)) { nna++; continue; }
        if (pos < 0 && !wneg) {
            wneg = 1;
            if (!quiet) fprintf(stderr, "physmerge: warning: negative position(s) found; "
                                        "block boundaries are clamped at 0.\n");
        }
        const char *id = (i_id >= 0 && fld[i_id]) ? fld[i_id] : "";
        nkept++;

        if (dosort) {
            int rank = -1;
            for (int k = 0; k < n_cn; k++) if (!strcmp(cnames[k], ch)) { rank = k; break; }
            if (rank < 0) { cnames = xrealloc(cnames, (size_t)(n_cn + 1) * sizeof(char *)); cnames[n_cn] = strdup(ch); rank = n_cn++; }
            if (n_rec == cap_rec) { cap_rec = cap_rec ? cap_rec * 2 : 65536; recs = xrealloc(recs, cap_rec * sizeof(Rec)); }
            recs[n_rec].chrom_rank = rank; recs[n_rec].pos = pos; recs[n_rec].val = val;
            recs[n_rec].idx = (unsigned long)n_rec;
            recs[n_rec].id_off = c.haveid ? aput(&arena, id) : (size_t)-1;
            recs[n_rec].line_off = annot ? aput(&arena, lcopy) : (size_t)-1;
            n_rec++;
            continue;
        }

        if (pset) {
            int same = c.havech ? (strcmp(sget(&pch), ch) == 0) : 1;
            if (same) {
                if (pos < ppos)
                    die("input is not position-sorted (chromosome %s: %.0f after %.0f).\n"
                        "       re-run with --sort, or sort the file first.", ch, pos, ppos);
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
        sset(&pch, ch); ppos = pos; pset = 1;

        push(&c, ch, pos, val, id, annot ? lcopy : "");
    }
    rd_close(&rd);

    if (dosort) {
        qsort(recs, n_rec, sizeof(Rec), rcmp);
        for (size_t k = 0; k < n_rec; k++)
            push(&c, cnames[recs[k].chrom_rank], recs[k].pos, recs[k].val,
                      aget(&arena, recs[k].id_off), aget(&arena, recs[k].line_off));
    }
    endch(&c);
    puthdr(&c);

    if (c.out != stdout) fclose(c.out);
    if (c.snpf) fclose(c.snpf);
    g_out = g_snp = NULL; g_snpf[0] = '\0';

    if (!quiet) {
        if (tfilt) fprintf(stderr, "physmerge: TEST filter: kept %lu of %lu rows where %s = '%s'.\n",
                                 nread - ntdrop, nread, tcol, tval);
        if (ncdrop) fprintf(stderr, "physmerge: %lu row(s) dropped by the chromosome filter.\n", ncdrop);
        if (nna) fprintf(stderr, "physmerge: %lu row(s) dropped (NA in position or value).\n", nna);
        fprintf(stderr, "physmerge: %lu SNPs -> %ld blocks (window=%.0f, sig_th=%g, reward=%s, reset_on=%s).\n",
                nkept, c.nblk, window, sig, rmax ? "max" : "min", rany ? "any" : "best");
    }
    return 0;
}
