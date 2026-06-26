 /*
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
/* math.h defines a macro isinf,
   this is omitted(!) if
   compiled with flag -ansi */
#include <math.h>
#include <float.h>
#include "itsolver.h"

#undef ITSOLVER_DEBUG

#define max(A,B) ((A) > (B) ? (A) : (B))

void convplot_steps(const char *filename,
  const char *title, const char *xlabel, const char *ylabel,
  const int k, const int *m, double **y, char **name)
{
  int i, j, ell, maxm = 0;
  int inc, incs, ifactor, istarty;
  int iminy, imaxy, ispany;
  double miny, maxy, xscale, yscale, color;
  FILE *out;
  time_t tt = time(NULL);
  struct tm* tmval = localtime(&tt);
  int year   = tmval->tm_year + 1900;
  int month  = tmval->tm_mon + 1;
  int day    = tmval->tm_mday;
  int hour   = tmval->tm_hour;
  int minute = tmval->tm_min;
  int second = tmval->tm_sec;

  /* open file */

  out = fopen(filename, "w");
  assert(out != NULL);

  /* print header, beware that % has to be %% */

  (void) fprintf(out,
    "%%!PS-Adobe-3.0 EPSF-3.0\n"
    "%%%%Creator: convplot_steps von Jens-Peter M. Zemke\n"
    "%%%%Title: Abbildung erzeugt von convplot_steps\n"
    "%%%%CreationDate: %d-%d-%d %d-%d-%d\n"
    "%%%%DocumentData: Clean7Bit\n"
    "%%%%Origin: 0 0\n"
    "%%%%BoundingBox: 0 0 600 400\n"
    "%%%%HiResBoundingBox: 0.0 0.0 600.0 400.0\n"
    "%%%%LanguageLevel: 2\n"
    "%%%%Pages: 1\n"
    "%%%%Page: 1 1\n",year,month,day,hour,minute,second);

  /* define handy commands */

  (void) fprintf(out,
    "/m {newpath moveto} bind def\n"
    "/s {setlinewidth} bind def\n"
    "/l {lineto} bind def\n"
    "/ls {lineto stroke} bind def\n"
    "/x {dup -0.25 m 0 ls} bind def\n"
    "/y {dup -0.25 exch m 0 exch ls} bind def\n"
    "/xt {dup -0.5 m 0 ls} bind def\n"
    "/yt {dup -0.5 exch m 0 exch ls} bind def\n"
    "/max {2 copy lt {exch} if pop} bind def\n");

  /* print scaling, translation and color information */

  (void) fprintf(out,
    "10 dup scale\n"
    "8 5 translate\n"
    "0 setgray\n");

  /* print title; xlabel; ylabel */

  (void) fprintf(out,
    "/Helvetica findfont 1.0 scalefont setfont\n"
    "(%s)\n"
    "dup 25 exch stringwidth pop 2 div sub 32.4 moveto show\n"
    "(%s)\n"
    "dup 25 exch stringwidth pop 2 div sub -3 moveto show\n"
    "(%s) dup 15 exch\n"
    "stringwidth pop 2 div sub -5 exch moveto\n"
    "90 rotate show -90 rotate\n", title, xlabel, ylabel);

  /* print axes */

  (void) fprintf(out,
    "0.08 s\n"
    "0 30 m\n"
    "0 dup l\n"
    "48 0 ls\n");
  /* "48 0 l 48 30 l 0 30 ls\n"); */

  /* compute size of axes */

  miny = DBL_MAX; /* largest positive double number */
  maxy = DBL_MIN*DBL_EPSILON; /* smallest positive double number */

  for (i=0; i<k; i++)
  {
    for (j=0; j<*m+1; j++)
    {
      if (!isinf(**y)) /* look for +/-inf */
      {
        if (**y > 0) /* only use positive vals */
        {
          if (miny > **y) miny = **y;
          if (maxy < **y) maxy = **y;
        }
      }
      (*y)++;
    }
    *y -= *m+1;
    if (maxm < *m) maxm = *m;
    m++;y++;
  }
  y -= k;
  m -= k;

#ifdef ITSOLVER_DEBUG
  printf("miny = %g, maxy = %g, maxm = %d\n", miny, maxy, maxm);
#endif

  /* set up semilog scale in y-direction */

  imaxy = ceil( log10(maxy));
  iminy = floor(log10(miny));
  ispany = imaxy-iminy;

  /* plot (almost) constant functions on top */

  if (ispany == 0)
  {
    iminy--;
    ispany++;
  }

#ifdef ITSOLVER_DEBUG
  printf("iminy = %d, imaxy = %d, ispany = %d\n",
    iminy, imaxy, ispany);
#endif

  /* compute scales for x- and y-ticks */

  xscale = 48.0/(maxm);
  yscale = 30.0/(imaxy-iminy);

  /* ticks and numbers on x-axis */

  ifactor = 1;
  inc = max(pow(10.0,round(log10(maxm))-1),1);
  incs = max(pow(10.0,round(log10(maxm))-2),1);

#ifdef ITSOLVER_DEBUG
  printf("x-axis: %d %d\n", inc, incs);
#endif

  if (maxm/inc > 18)
  {
    ifactor = 5;
    incs = inc;
  }
  else if (maxm/inc > 10)
  {
    ifactor = 2;
  }

  (void) fprintf(out,"0.05 s\n");

  for (i=0; i<=maxm; i += ifactor*inc)
  {
    (void) fprintf(out,"%.3f xt ", i*xscale);
    for (j=incs; j<ifactor*inc; j += incs)
    {
      if (i+j<=maxm)
      {
        (void) fprintf(out,"%.3f x ", (i+j)*xscale);
      }
    }
    (void) fprintf(out, "\n"
      "/Helvetica findfont 1.0 scalefont setfont\n"
      "(%d) dup %.3f exch stringwidth pop 2 div "
      "sub -1.4 moveto show\n", i, i*xscale);
  }

  /* ticks and numbers on y-axis */

  ifactor = 1;
  inc = max(pow(10.0,floor(log10(ispany/2))),1);
  incs = max(pow(10.0,floor(log10(ispany/2))-1),1);

#ifdef ITSOLVER_DEBUG
  printf("y-axis: %d %d\n", inc, incs);
#endif

  if (ispany/inc > 20)
  {
    ifactor = 5;
    incs = inc;
  }
  if (ispany/inc > 10)
  {
    ifactor = 2;
  }

  /* length of strings in exponents in y-axis */
  
  (void) fprintf(out,"/yl 0 def\n"
    "/Helvetica findfont 0.8 scalefont setfont\n");
  (void) fprintf(out,
    "/yl (%d) stringwidth pop yl max def\n",iminy);
  (void) fprintf(out,
    "/yl (%d) stringwidth pop yl max def\n",imaxy);

  (void) fprintf(out,"0.05 s\n");

  /* print fine y-ticks below first mark (if any) */

  istarty = -iminy%inc;

  for (ell = istarty-incs; ell >= 0; ell -= incs)
  {
    (void) fprintf(out,"%.3f y ", ell*yscale);
  }
  if (istarty-incs > 0)
  {
    (void) fprintf(out,"\n");
  }

  /* print major y-ticks with fine y-ticks */

  for (i=istarty; i<=ispany; i += ifactor*inc)
  {
    (void) fprintf(out,"%.3f yt\n", i*yscale);
    for (j=incs; j<ifactor*inc; j += incs)
    {
      if (i+j<=ispany)
      {
        (void) fprintf(out,"%.3f y\n", (i+j)*yscale);
      }
    }
    (void) fprintf(out,
      "/Helvetica findfont 1.0 scalefont setfont\n"
      "-1.95 yl sub %.3f  moveto (10) show\n", i*yscale-.5);
    (void) fprintf(out,
      "/Helvetica findfont 0.8 scalefont setfont\n"
      "-.7 yl sub %.3f moveto (%d) show\n", i*yscale, i+iminy);
  }

  /* compute maximal length of strings in legend */
  
  (void) fprintf(out,"/nl 0 def\n"
    "/Helvetica findfont 1.0 scalefont setfont\n");
  for (i=0; i<k; i++)
  {
    (void) fprintf(out,
       "/nl (%s) stringwidth pop nl max def\n",name[i]);
  }

  /* select linewidth */

  (void) fprintf(out,"0.08 s\n");

  for (i=0; i<k; i++)
  {

    /* select color, move to initial poition */

    color = (3.0*i)/k;

    switch((int) floor(color))
    {
      case 0:
        (void) fprintf(out,
          "%f 0 %f setrgbcolor\n",color,1.0-color);
        break;
      case 1:
        (void) fprintf(out,
	  "%f %f 0 setrgbcolor\n",2.0-color,color-1.0);
        break;
      case 2:
        (void) fprintf(out,
	  "0 %f %f setrgbcolor\n",3.0-color,color-2.0);
        break;
    }

    /* print convergence curve */

    j = 0;
    while ( (y[i][j] <= 0.0) || isinf(y[i][j])) j++;
    (void) fprintf(out,"%6.3f %6.3f m\n",
      j*xscale, (log10(y[i][j])-iminy)*yscale);

    for (ell=j; ell<m[i]+1; ell++)
    {
      if ( (y[i][ell] > 0.0) && (!isinf(y[i][ell])) )
      {
        (void) fprintf(out,"%6.3f %6.3f l\n",
          ell*xscale, (log10(y[i][ell])-iminy)*yscale);
      }
    }
    (void) fprintf(out,"stroke\n");

    /* generate line and name in legend */

    (void) fprintf(out,
      "48 nl 2.4 add sub %.3f moveto\n"
      "1.5 0 rlineto stroke\n"
      "48 nl 0.4 add sub %.3f moveto\n"
      "/Helvetica findfont 1.0 scalefont setfont\n"
      "0 setgray (%s) show\n",
      29.0-i+.35, 29.0-i, name[i]);
  }

  /* print border around legend */

  (void) fprintf(out,
    "48 30 m\n"
    "nl 2.9 add neg 0 rlineto\n"
    "0 %.3f rlineto\n"
    "nl 2.9 add 0 rlineto\n"
    "closepath\n"
    "stroke\n",-k*1.0-.3); /* 2*.65 = 1.3 */

  /* print footer */

  (void) fprintf(out,"showpage\n");
  (void) fclose(out);

} /* convplot_steps */
