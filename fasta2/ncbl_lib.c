/*	ncbl_lib.c	functions to read ncbi-blast format files from
			setdb (blastp 1.3.2) format files

		copyright (c) 1992 William R. Pearson
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef UNIX
#define RBSTR "r"
#else
#define RBSTR "rb"
#endif

#define XTERNAL
#include "uascii.gbl"
#include "upam.gbl"
#include "ncbl_head.h"

extern int ldnaseq;
void src_ulong_read();
void src_char_read();

/* nt_btoa maps  from blast 2bit  format to ascii  characters */
static char nt_btoa[5] = {"ACGT"};

static char aa_btoa[27]= {"-ARNDCQEGHILKMFPSTWYVBZX*"};
static int aa_btof[32];	/* maps to fasta alphabet */

static FILE *tfile=NULL,	/* table of offsets, also DB info */
	    *hfile=NULL,	/* description lines */
	    *sfile=NULL;	/* binary sequence data */

static unsigned long lib_cnt, max_cnt, totlen, mxlen, dbline_len;
static unsigned long *seq_beg, *hdr_beg;
static unsigned char *ambiguity_ray;
static long seq_format, dbtype, dbformat;
static char dline[512];

#define NCBIBL13 11

extern int (*getlib)();
int ncbl_getliba(), ncbl_getlibn();

int ncbl_openlib(name)
     char *name;
{
  char hname[256];
  char sname[256];
  char tname[256];
  long title_len;
  char *title_str;
  int rdtmp;
  int i;
  unsigned long line_len, c_len, clean_count;

  if (!ldnaseq) {
    newname(tname,name,AA_TABLE_EXT,(int)sizeof(tname));
    if ((tfile = fopen(tname,RBSTR))==NULL) {
      fprintf(stderr," cannot open %s (%s.%s) table file\n",
	      name,tname,NT_TABLE_EXT);
      return (-1);
    }
    seq_format = AAFORMAT;
  }
  else {
    newname(tname,name,NT_TABLE_EXT,(int)sizeof(tname));
    if ((tfile = fopen(tname,RBSTR))==NULL) {
      fprintf(stderr," cannot open %s (%s.%s) table file\n",
	      name,tname,NT_TABLE_EXT);
      return (-1);
    }
    seq_format = NTFORMAT;
  }
	
  src_ulong_read(tfile,&dbtype);
  src_ulong_read(tfile,&dbformat);

  if (seq_format == AAFORMAT && (dbformat != seq_format || dbtype !=
				 DB_TYPE_PRO)) {
    fprintf(stderr,"error - %s wrong type (%ld/%d) or format (%ld/%ld)\n",
	    tname,dbtype,DB_TYPE_PRO,dbformat,seq_format);
    return (-1);
  }
  else if (seq_format == NTFORMAT && (dbformat != seq_format || dbtype !=
				 DB_TYPE_NUC)) {
    fprintf(stderr,"error - %s wrong type (%ld/%d) or format (%ld/%ld)\n",
	    tname,dbtype,DB_TYPE_NUC,dbformat,seq_format);
    return (-1);
  }

  if (seq_format == AAFORMAT) {
    newname(hname,name,AA_HEADER_EXT,(int)sizeof(hname));
    if ((hfile = fopen(hname,RBSTR))==NULL) {
      fprintf(stderr," cannot open %s header file\n",hname);
      return (-1);
    }
    newname(sname,name,AA_SEARCHSEQ_EXT,(int)sizeof(sname));
    if ((sfile = fopen(sname,RBSTR))==NULL) {
      fprintf(stderr," cannot open %s sequence file\n",sname);
      return (-1);
    }
    getlib = ncbl_getliba;
  }
  else {
    newname(hname,name,NT_HEADER_EXT,(int)sizeof(hname));
    if ((hfile = fopen(hname,RBSTR))==NULL) {
      fprintf(stderr," cannot open %s header file\n",hname);
      return (-1);
    }
    newname(sname,name,NT_SEARCHSEQ_EXT,(int)sizeof(sname));
    if ((sfile = fopen(sname,RBSTR))==NULL) {
      fprintf(stderr," cannot open %s sequence file\n",sname);
      return (-1);
    }
    getlib = ncbl_getlibn;
  }

/* all files should be open */

  src_ulong_read(tfile,&title_len);
  rdtmp = title_len + ((title_len%4 !=0 ) ? 4-(title_len%4) : 0);
  if ((title_str = calloc((size_t)rdtmp,sizeof(char)))==NULL) {
    fprintf(stderr," cannot allocate title string (%d)\n",rdtmp);
    return(-1);
  }
  fread(title_str,(size_t)1,(size_t)rdtmp,tfile);

  lib_cnt = 0;
  if (seq_format == AAFORMAT) {
    src_ulong_read(tfile,&max_cnt);
    src_ulong_read(tfile,&totlen);
    src_ulong_read(tfile,&mxlen);

    /* fprintf(stderr," max_cnt: %d, totlen: %d\n",max_cnt,totlen); */

    if ((seq_beg=(unsigned long *)calloc((size_t)max_cnt+1,sizeof(long)))==NULL) {
      fprintf(stderr," cannot allocate sequence pointers\n");
      return -1;
    }
    if ((hdr_beg=(unsigned long *)calloc((size_t)max_cnt+1,sizeof(long)))==NULL) {
      fprintf(stderr," cannot allocate header pointers\n");
      return -1;
    }
    for (i=0; i<max_cnt+1; i++) src_ulong_read(tfile,&seq_beg[i]);
    for (i=0; i<max_cnt+1; i++) src_ulong_read(tfile,&hdr_beg[i]);

    for (i=0; i<sizeof(aa_btoa); i++) {
      if ((rdtmp=aascii[aa_btoa[i]])<NA) aa_btof[i]=rdtmp;
      else aa_btof[i]=aascii['X'];
    }
  }
  else if (seq_format == NTFORMAT) {
    src_ulong_read(tfile,&dbline_len);	/* length of uncompress DB lines */
    src_ulong_read(tfile,&max_cnt);	/* number of entries */
    src_ulong_read(tfile,&mxlen);	/* maximum length sequence */
    src_ulong_read(tfile,&totlen);	/* total count */
    src_ulong_read(tfile,&c_len);	/* compressed db length */
    src_ulong_read(tfile,&clean_count);	/* count of nt's cleaned */

    fseek(tfile,(size_t)((clean_count)*4),1);
    					 /* seek over clean_count */
    if ((seq_beg=(unsigned long *)calloc((size_t)max_cnt+1,sizeof(long)))==NULL) {
      fprintf(stderr," cannot allocate sequence pointers\n");
      return -1;
    }
    if ((hdr_beg=(unsigned long *)calloc((size_t)max_cnt+1,sizeof(long)))==NULL) {
      fprintf(stderr," cannot allocate header pointers\n");
      return -1;
    }
    if ((ambiguity_ray=
	 (unsigned char *)calloc((size_t)max_cnt/CHAR_BIT+1,sizeof(char)))==NULL) {
      fprintf(stderr," cannot allocate ambiguity_ray\n");
      return -1;
    }

    for (i=0; i<max_cnt+1; i++) src_ulong_read(tfile,&seq_beg[i]);
    fseek(tfile,(size_t)((max_cnt+1)*4),1);
    					 /* seek over seq_beg */
    for (i=0; i<max_cnt+1; i++) src_ulong_read(tfile,&hdr_beg[i]);
    for (i=0; i<max_cnt/CHAR_BIT+1; i++)
      src_char_read(tfile,&ambiguity_ray[i]);
  }
  return 1;
}

void ncbl_closelib()
{
  if (tfile !=NULL ) {fclose(tfile); tfile=NULL;}
  if (hfile !=NULL ) {fclose(hfile); hfile=NULL;}
  if (sfile !=NULL ) {fclose(sfile); sfile=NULL;}
}

int ncbl_getliba(seq,maxs,libstr,libpos,lcont)
     char *seq;
     int maxs;
     char *libstr;
     long *libpos;
     int *lcont;
{
  register char *sptr;
  long seqcnt;
  long tmp;
  char ch;
  static long seq_len;
  
  *libpos = lib_cnt;
  if (*lcont==0) {
    if (lib_cnt >= max_cnt) return -1;
    seq_len = seq_beg[lib_cnt+1] - seq_beg[lib_cnt] -1;
    tmp=fgetc(sfile);	/* skip the null byte */
    if (tmp!=NULLB) fprintf(stderr," phase error: %ld:%ld found\n",lib_cnt,tmp);
    libstr[0]='\0';
    }
  
  if (seq_len < maxs) {
    if ((tmp=fread(seq,(size_t)1,(size_t)seq_len,sfile))!=(size_t)seq_len) {
      fprintf(stderr," could not read sequence record: %ld %ld != %ld\n",
	      *libpos,tmp,seq_len);
      goto error; 
    }
    if (aa_btoa[seq[seq_len-1]]=='*') seqcnt = seq_len-1;
    else seqcnt=seq_len;
    lib_cnt++;
    *lcont = 0;
  }
  else {
    if (fread(seq,(size_t)1,(size_t)(maxs-1),sfile)!=(size_t)(maxs-1)) {
      fprintf(stderr," could not read sequence record: %ld %ld\n",
	      *libpos,seq_len);
      goto error;
    }
    (*lcont)++;
    seqcnt = maxs-1;
    seq_len -= seqcnt;
  }
  sptr = seq+seqcnt;

  while (--sptr >= seq) *sptr = aa_btof[*sptr];
  
  seq[seqcnt]= EOSEQ;
  return (seqcnt);
  
error:	fprintf(stderr," error reading %ld at %ld\n",libstr,*libpos);
  fflush(stderr);
  return (-1);
}

ncbl_getlibn(seq,maxs,libstr,libpos,lcont)
     char *seq;
     int maxs;
     char *libstr;
     long *libpos;
     int *lcont;
{
  register char *sptr, *tptr, stmp;
  long seqcnt;
  long tmp;
  char ch;
  static long seq_len;
  static int c_len,c_pad;
  
  *libpos = lib_cnt;
  if (*lcont==0) {
    if (lib_cnt >= max_cnt) return -1;
    c_len = seq_beg[lib_cnt+1]/(CHAR_BIT/NBPN)
		- seq_beg[lib_cnt]/(CHAR_BIT/NBPN);
    c_len -= NSENTINELS;

    seq_len = c_len*(CHAR_BIT/NBPN);
    c_pad = seq_beg[lib_cnt] & ((CHAR_BIT/NBPN)-1);
    if (c_pad != 0) seq_len -= ((CHAR_BIT/NBPN) - c_pad);

    tmp=fgetc(sfile);	/* skip the null byte */
    if (tmp!=NT_MAGIC_BYTE) {
      fprintf(stderr," phase error: %ld:%ld (%ld/%d) found\n",
	      lib_cnt,seq_len,tmp,NT_MAGIC_BYTE);
      goto error;
    }
    libstr[0]='\0';
  }

  if (seq_len < maxs-3) {
    seqcnt=(seq_len+3)/4;
    if (seqcnt==0) seqcnt++;
    if ((tmp=fread(seq,(size_t)1,(size_t)seqcnt,sfile))
	!=(size_t)seqcnt) {
      fprintf(stderr," could not read sequence record: %s %ld %ld != %ld: %d\n",
	      libstr,*libpos,tmp,seqcnt,*seq);
      goto error; 
    }
    tmp=fgetc(sfile);	/* skip the null byte */
    if (tmp!=(unsigned char)NT_MAGIC_BYTE) {
      fprintf(stderr," phase2 error: %ld:%ld (%ld/%d) next ",
	      lib_cnt,seqcnt,tmp,NT_MAGIC_BYTE);
      
      goto error;
    }
    *lcont = 0;
    lib_cnt++;
  }
  else {
    seqcnt = ((maxs+3)/4)-1;
    if (fread(seq,(size_t)1,(size_t)(seqcnt),sfile)!=(size_t)(seqcnt)) {
      fprintf(stderr," could not read sequence record: %s %ld %ld\n",
	      libstr,*libpos,seqcnt);
      goto error;
    }
    (*lcont)++;
  }
  
  /* point to the last packed byte and to the end of the array
     seqcnt is the exact number of bytes read
     tptr points to the destination, use multiple of 4 to simplify math
     sptr points to the source, note that the last byte will be read 4 cycles
     before it is written
     */
  
  sptr = seq + seqcnt;
  tptr = seq + 4*seqcnt;
  while (sptr>seq) {
    stmp = *--sptr;
    *--tptr = stmp&3;
    *--tptr = (stmp >>= 2)&3;
    *--tptr = (stmp >>= 2)&3;
    *--tptr = (stmp >>= 2)&3;
  }
  /*
    for (sptr=seq; sptr < seq+seq_len; sptr++) {
    printf("%c",nt[*sptr]);
    if ((int)(sptr-seq) % 60 == 59) printf("\n");
    }
    printf("\n");
    */
  if (seqcnt*4 >= seq_len) {	/* there was enough room */
    seq[seq_len]= EOSEQ;
    /* printf("%d\n",seq_len); */
    return seq_len;
  }
  else {				/* not enough room */
    seq[seqcnt*4]=EOSEQ;
    seq_len -= 4*seqcnt;
    return (4*seqcnt);
  }
  
error:	fprintf(stderr," error reading %ld at %ld\n",libstr,*libpos);
  fflush(stderr);
  return (-1);
}

ncbl_ranlib(str,cnt,libpos)
	char *str; int cnt;
	long libpos;
{
  char hline[256], *bp, *bp0;
  int llen;
  long spos;

  lib_cnt = libpos;
  llen = hdr_beg[lib_cnt+1]-hdr_beg[lib_cnt];
  if (llen > sizeof(hline)) llen = sizeof(hline);
  fseek(hfile,hdr_beg[lib_cnt]+1,0);

  fread(hline,(size_t)1,(size_t)(llen-1),hfile);
  hline[llen-1]='\0';

  if (hline[9]=='|' || hline[10]=='|') {
    bp0 = strchr(hline+3,'|');
    if ((bp=strchr(bp0+1,' '))!=NULL) *bp='\0';
    if (dbformat == NTFORMAT && 
	(ambiguity_ray[lib_cnt/CHAR_BIT]&(1<<lib_cnt%CHAR_BIT))) {
      sprintf(str,"*%-9s ",bp0+1);
    }
    else sprintf(str,"%-10s ",bp0+1);
    strncat(str+11,bp+1,cnt-strlen(str));
  }
  else {
    if (dbformat == NTFORMAT && 
	(ambiguity_ray[lib_cnt/CHAR_BIT]&(1<<lib_cnt%CHAR_BIT))) {
      str[0]='*'; 
      strncpy(str+1,hline,cnt-1);
    }
    else strncpy(str,hline,cnt);
  }
  str[cnt-1]='\0';

  if (dbformat == AAFORMAT)
    fseek(sfile,seq_beg[lib_cnt]-1,0);
  else {
    spos = (seq_beg[lib_cnt])/(CHAR_BIT/NBPN);
    fseek(sfile,spos-1,0);
  }
}

void src_ulong_read(fd, val)
     FILE *fd;
     unsigned long *val;
{
#ifdef BIG_ENDIAN
  fread((char *)val,(size_t)4,(size_t)1,fd);
#else
  unsigned char b[4];

  fread((char *)&b[0],(size_t)1,(size_t)4,fd);
  *val = 0;
  *val = (unsigned long)((unsigned long)((unsigned long)(b[0]<<8) +
	 (unsigned long)b[1]<<8) + (unsigned long)b[2]<<8)+(unsigned long)b[3];
#endif
}

void src_long_read(fd,val)
     FILE *fd;
     long *val;
{
#ifdef BIG_ENDIAN
  fread((char *)val,(size_t)4,(size_t)1,fd);
#else
  unsigned char b[4];

  fread((char *)&b[0],(size_t)1,(size_t)4,fd);
  *val = 0;
  *val = (long)((long)((long)(b[0]<<8)+(long)b[1]<<8)+(long)b[2]<<8)
	  +(long)b[3];
#endif
}

void src_char_read(fd, val)
     FILE *fd;
     char *val;
{
  fread(val,(size_t)1,(size_t)1,fd);
}

void src_fstr_read(fd, val, slen)
     FILE *fd;
     char *val;
     long slen;
{
  fread(val,(size_t)slen,(size_t)1,fd);
}

newname(nname,oname,suff,maxn)
	char *nname, *oname, *suff;
{
	char *tptr;

	if (oname[0]=='@') strncpy(nname,&oname[1],maxn);
	else strncpy(nname,oname,maxn);
	for (tptr=nname; *tptr=='.' && *tptr; tptr++);
	for (; *tptr!='.'&& *tptr; tptr++); /* get to '.' or EOS */
	*tptr++='.'; *tptr='\0';
	strncat(nname,suff,maxn);
	}

