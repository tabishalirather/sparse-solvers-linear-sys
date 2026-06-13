 /*
 *  \author Michael Wende
 *  \author Jens-Peter Zemke
 *  \author Sabine Le Borne
 */

#include "itsolver.h"

struct conv_plot_line_s;
typedef struct conv_plot_line_s conv_plot_line_t;

struct conv_plot_line_s
{
  /* Number of steps performed */
  int steps;

  /* Vector of (steps + 1) residual values to be plotted. */
  double *y;

  /* Legend entry for this line of the plot.
   * Each entry may consist of at most LEGEND_ENTRY_LENGTH characters. */
  char *name;

  /* Next and previous lines in circular and doubly linked list. */
  conv_plot_line_t *next;
  conv_plot_line_t *prev;
};

/* Maximum string length for the plot title. */
#define TITLE_LENGTH (100)

/* Maximum string length for legend entries. */
#define LEGEND_ENTRY_LENGTH (40)

/* Initialize dummy head corresponding to an empty list.
 */
static void conv_plot_line_init(conv_plot_line_t *head);

/* Create a new line and append it to the circular list with the given empty
 * head.
 * Make copies of the given residual values and legend name before storing
 * them in the line structure so that these may be freed regardless of their
 * origin when the line is destroyed.
 */
static void conv_plot_line_append(conv_plot_line_t *head, int steps, const double *y, char *name);

/* Destroy a line and remove it from the list in which it resides.
 */
static void conv_plot_line_delete(conv_plot_line_t *conv_plot_line);

/*! \brief Convergence plot structure. */
struct conv_plot_s
{
  /*! \brief Number of lines in this plot. */
  int num;

  /*! \brief Plot title, copy of the string passed in. */
  char *title;

  /*! \brief Horizontal axis descriptions. */
  char *xlabel;

  /*! \brief Vertical axis descriptions. */
  char *ylabel;

  /*! \brief Empty head of a circular linked list where each entry corresponds to a
   *         line of the plot. */
  conv_plot_line_t head;
};

static void conv_plot_line_init(conv_plot_line_t *head)
{
  assert(head);
  head->steps = -1;
  head->y = 0;
  head->name = 0;
  head->next = head->prev = head;
}

static void conv_plot_line_append(conv_plot_line_t *head, int steps, const double *y, char *name)
{
  assert(head);
  assert(steps >= 0);
  assert(y);
  conv_plot_line_t *line = malloc(sizeof(conv_plot_line_t));
  line->steps = steps;
  line->y = malloc((steps + 1) * sizeof(double));
  vec_copy(steps + 1, y, line->y);
  line->name = name ? strdup(name) : 0;
  line->next = head;
  line->prev = head->prev;
  head->prev->next = line;
  head->prev = line;

}

static void conv_plot_line_delete(conv_plot_line_t *line)
{
  /* Remove line from its list. */
  line->prev->next = line->next;
  line->next->prev = line->prev;

  /* Delete structure and its members. */
  free(line->y);
  free(line->name);
  memset(line, 0, sizeof(conv_plot_line_t));
  free(line);
}

conv_plot_t *conv_plot_new()
{
  conv_plot_t *conv_plot = malloc(sizeof(conv_plot_t));
  conv_plot->num = 0;
  conv_plot->title = 0;
  conv_plot->xlabel = 0;
  conv_plot->ylabel = 0;
  conv_plot_line_init(&conv_plot->head);
  return conv_plot;
}

void conv_plot_title(conv_plot_t *conv_plot, const char *title)
{
  assert(conv_plot);
  assert(title);
  assert(!conv_plot->title);
  conv_plot->title = strndup(title, TITLE_LENGTH);
}

void conv_plot_titlef(conv_plot_t *conv_plot, const char *title, ...)
{
  /* Assert that no title has been set yet. */
  assert(conv_plot);
  assert(title);
  assert(!conv_plot->title);

  /* Set new title. */
  va_list ap;
  va_start(ap, title);
  conv_plot->title = malloc(TITLE_LENGTH);
  vsnprintf(conv_plot->title, TITLE_LENGTH, title, ap);
  va_end(ap);
}

void conv_plot_xlabel(conv_plot_t *conv_plot, const char *xlabel)
{
  assert(conv_plot);
  assert(xlabel);
  assert(!conv_plot->xlabel);
  conv_plot->xlabel = strdup(xlabel);
}

void conv_plot_ylabel(conv_plot_t *conv_plot, const char *ylabel)
{
  assert(conv_plot);
  assert(ylabel);
  assert(!conv_plot->ylabel);
  conv_plot->ylabel = strdup(ylabel);
}

void conv_plot_line(conv_plot_t *conv_plot, int steps, const double *y)
{
  assert(conv_plot);
  assert(steps >= 0);
  assert(y);
  ++conv_plot->num;
  conv_plot_line_append(&conv_plot->head, steps, y, 0);
}

void conv_plot_name(conv_plot_t *conv_plot, const char *name)
{
  assert(conv_plot);
  assert(name);
  conv_plot_line_t *last = conv_plot->head.prev;
  assert(last != &conv_plot->head);
  assert(!last->name);

  /* Allocate a legend entry string of a limited length.
   * Any characters exceeding this length will be dropped. */
  last->name = strndup(name, LEGEND_ENTRY_LENGTH);
}

void conv_plot_namef(conv_plot_t *conv_plot, const char *name, ...)
{
  assert(conv_plot);
  assert(name);
  conv_plot_line_t *last = conv_plot->head.prev;
  assert(last != &conv_plot->head);
  assert(!last->name);

  /* Allocate a legend entry string of a limited length.
   * Any characters exceeding this length will be dropped. */
  last->name = malloc(LEGEND_ENTRY_LENGTH);

  /* Pass va_list to vsnprintf. */
  va_list ap;
  va_start(ap, name);
  vsnprintf(last->name, LEGEND_ENTRY_LENGTH, name, ap);
  va_end(ap);
}

void conv_plot_save(conv_plot_t *conv_plot, const char *file_name)
{
  assert(conv_plot);
  assert(file_name);

  /* Axis limits cannot be computed when the plot is empty. */
  assert(conv_plot->num > 0);

  /* Add .eps suffix to file name. */
  char *file_name_ps = malloc(strlen(file_name) + 4 + 1);
  sprintf(file_name_ps, "%s.eps", file_name);

  /* Collect vectors. */
  int *m = malloc(conv_plot->num * sizeof(int));
  double **y = malloc(conv_plot->num * sizeof(double *));
  char **name = malloc(conv_plot->num * sizeof(char *));
  conv_plot_line_t *line = conv_plot->head.next;
  int i = 0;
  while (line != &conv_plot->head) {
    assert(i < conv_plot->num);
    m[i] = line->steps;
    y[i] = line->y;
    name[i] = line->name ? line->name : "";
    line = line->next;
    ++i;
  }
  assert(i == conv_plot->num);

  /* Call Jens. */
  convplot_steps(file_name_ps,
                 conv_plot->title ? conv_plot->title : "",
                 conv_plot->xlabel ? conv_plot->xlabel : "Step",
                 conv_plot->ylabel ? conv_plot->ylabel : "Residual",
     conv_plot->num,
     m, y, name);

  /* Clean up. */
  free(file_name_ps);
  free(m);
  free(y);
  free(name);
}

void conv_plot_delete(conv_plot_t *conv_plot)
{
  assert(conv_plot);
  free(conv_plot->title);
  free(conv_plot->xlabel);
  free(conv_plot->ylabel);
  while (conv_plot->head.next != &conv_plot->head) {
    conv_plot_line_delete(conv_plot->head.next);
  }
  memset(conv_plot, 0, sizeof(conv_plot_t));
  free(conv_plot);
}
