/* 
   qd - qd is a simple, pedagogical implementation of TDDFT for 2D systems.

   Copyright (C) 2010 The 2010 Benasque TDDFT School
*/

#include <stdio.h>
#include <sys/types.h>
#include <argp.h>
#include "system.h"

#define EXIT_FAILURE 1

#define PACKAGE "QDOT"
#define VERSION "0.1"

#if ENABLE_NLS
# include <libintl.h>
# define _(Text) gettext (Text)
#else
# define textdomain(Domain)
# define _(Text) Text
#endif
#define N_(Text) Text

char *xmalloc ();
char *xrealloc ();
char *xstrdup ();

static error_t parse_opt (int key, char *arg, struct argp_state *state);
static void show_version (FILE *stream, struct argp_state *state);

/* argp option keys */
enum {DUMMY_KEY=129
};

/* Option flags and variables.  These are initialized in parse_opt.  */


static struct argp_option options[] =
{
  { "coefficients", 'c', 0,      0, "Generates coefficients for the discretization"},
  { "test_hartree", 'h', 0,      0, "Tests the Poisson solver"},
  { "test_laplacian", 'l', 0,      0, "Tests the Laplacian"},
  { "test_exponential", 'e', 0,      0, "Tests the exponential"},
  { "gs", 'g', 0,      0, "Performs a ground state calculation"},
  { "td", 't', 0,      0, "Performs a time-dependent calculation"},
  { "strength_function", 's', 0,      0, "Calculates the strength function"},
  { "excitations", 'x', 0,      0, "Performs a LR-TDDFT calculation"},
  { NULL, 0, NULL, 0, NULL, 0 }
};

/* The argp functions examine these global variables.  */
const char *argp_program_bug_address = "<alberto@physik.fu-berlin.de>";
void (*argp_program_version_hook) (FILE *, struct argp_state *) = show_version;

/* Used by `main' to communicate with `parse_opt'. */
struct arguments
{
  /*char *args[2]; */               /* ARG1 & ARG2 */
  /*int coefficients;*/
  int mode;
  /*char *output_file;*/
};



static struct argp argp =
{
  //  options, parse_opt, N_("[FILE...]"),
  options, parse_opt, NULL,
  N_("qd is a simple, pedagogical implementation of TDDFT for 2D systems."),
  NULL, NULL, NULL
};



/* Fortran wrappers */
#define FORTRANMAIN_FC FC_FUNC (fortranmain, FORTRANMAIN)
#ifdef __cplusplus
extern "C"  /* prevent C++ name mangling */
#endif
void FORTRANMAIN_FC(int *mode);


int
main (int argc, char **argv)
{
  textdomain(PACKAGE);
  struct arguments arguments;

  /* Default values. */
  arguments.mode = 0;
  /*arguments.verbose = 0;*/
  /*arguments.output_file = "-";*/
  argp_parse(&argp, argc, argv, 0, NULL, &arguments);

  show_version(stdout, NULL);

  /* This does the work */
  FORTRANMAIN_FC(&arguments.mode);

  exit (0);
}

/* Parse a single option.  */
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  /* Get the INPUT argument from `argp_parse', which we
     know is a pointer to our arguments structure. */
  struct arguments *arguments = state->input;

  switch (key)
    {
    case 'c':
      arguments->mode = 1;
      break;
    case 'h':
      arguments->mode = 2;
      break;
    case 'l':
      arguments->mode = 3;
      break;
    case 'e':
      arguments->mode = 4;
      break;
    case 'g':
      arguments->mode = 5;
      break;
    case 't':
      arguments->mode = 6;
      break;
    case 's':
      arguments->mode = 7;
      break;
    case 'x':
      arguments->mode = 8;
      break;

    case ARGP_KEY_INIT:
      /* Set up default values.  */
      break;


    case ARGP_KEY_ARG:		/* [FILE]... */
      /* TODO: Do something with ARG, or remove this case and make
         main give argp_parse a non-NULL fifth argument.  */
      break;

    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

/* Show the version number and copyright information.  */
static void
show_version (FILE *stream, struct argp_state *state)
{
  (void) state;
  /* Print in small parts whose localizations can hopefully be copied
     from other programs.  */
  fputs(PACKAGE" "VERSION"\n", stream);
  fprintf(stream, _("Written by %s.\n\n"), "The 2010 Benasque TDDFT School");
  fprintf(stream, _("Copyright (C) %s %s\n"), "2010", "The 2010 Benasque TDDFT School");
  fputs(_("\
This program is free software; you may redistribute it under the terms of\n\
the GNU General Public License.  This program has absolutely no warranty.\n"),
	stream);
}
