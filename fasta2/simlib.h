#ifndef SIMLIB_H
#define SIMLIB_H

typedef int bool;
typedef unsigned char uchar;

/* for reading sequence files : */
int get_seq(int seqnbr, char *file_name, /*@out@*/ uchar ** seqptr);
bool is_DNA(uchar *s, int len);

/* for getting command-line parameters */
void ckargs(char *options, int argc, char **argv, int non_options);
bool get_argval(char c, int *val_ptr);

/* for DNA substitution scores and the statistics-based threshold */
void DNA_scores( int *O_ptr, int *E_ptr, int **ss);
int DNA_thresh();

/* for protein substitution scores and the statistics-based threshold */
void protein_scores(int *O_ptr, int *E_ptr, int **ss);
int protein_thresh();

/* for printing alignments */
void print_align_header();
void print_align(int score, int beg1, int end1, int beg2, int end2,
	int *script);
void print_block(int score, int beg1, int end1, int beg2, int end2, int f);

/* utility routines */
#ifdef __GNUC__
#define NORETURN __attribute__((__noreturn__))
#else
#define NORETURN /* */
#endif

/*@exits@*/ void fatal(char *msg) NORETURN;
/*@exits@*/ void fatalf(char *msg, char *val) NORETURN;
FILE *ckopen(char *name, char *mode);
void *ckalloc(size_t amount);
bool strsame(char *s, char *t);
char *strsave(char *s);

#endif /* SIMLIB_H */
