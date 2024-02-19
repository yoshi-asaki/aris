#include <stdlib.h>
#include <stdio.h>


int  main(int argc, char **argv)
{
  static char  string[200];
  static int   i, j, k;
  static FILE  *ifp, *ofp;

  if (argc != 3) {
    printf("Usage: tout INFILE OUTFILE\n");
    exit (0);
  }

  ifp = fopen(argv[1], "r");
  ofp = fopen(argv[2], "w");

  while (1) {
    if (fgets(string, sizeof(string), ifp) == NULL) {
      break;
    }

    for (i=0; i<strlen(string); i++) {
      if (string[i] != ' ' && string[i] != ',' && string[i+1] == ' ') {
        string[i+1] = ',';
      } else if (string[i] == ',' && string[i+1] == '\n') {
        string[i+2] ='\0';
        break;
      } else if (string[i] != ' ' && string[i+1] == '\n') {
        string[i+1] = ',';
        string[i+2] = '\n';
        string[i+3] ='\0';
        break;
      } else if (string[i] == '\n') {
        string[i+1] ='\0';
        break;
      }
    }
    fprintf(ofp, "%s", string);

  }

  fclose (ifp);
  fclose (ofp);
}
