#include "pgm.h"
#include <assert.h>
#include <stdio.h>


int main (int ac, char **av) {  /* av[1] contient le nom de l'image, av[2] le nom du resultat . */
  int nl,nc, oldnl,oldnc;
  unsigned char** im1=NULL;
  double **final=NULL;
  if (ac < 2) {printf("Usage : %s entree sortie \n",av[0]); exit(1); }
  /* Lecture d'une image pgm dont le nom est passe sur la ligne de commande */
  im1=lectureimagepgm(av[1],&nl,&nc);
  if (im1==NULL)  { puts("Lecture image impossible"); exit(1); }
  /* Debruitage de l'image*/
  oldnl=nl;
  oldnc=nc;
  final=ContourDetection(im1,&nl,&nc);
  /* Sauvegarde dans un fichier dont le nom est passe sur la ligne de commande */
  ecritureimagepgm(av[2],crop(imdouble2uchar(final,nl,nc),0,0,oldnl,oldnc),oldnl,oldnc);
  return 0;
}
