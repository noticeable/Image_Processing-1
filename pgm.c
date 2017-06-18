#include "pgm.h"

static int match_key(int fd, char *key) { char buf[80];
  read(fd, buf, strlen(key));
  if( strncmp(buf, key, strlen(key)) != 0 ) return FALSE;
  else return TRUE;
}

static void skip_comment(int fd, char code, char *c) {
  while( *c == code ) {
      while( (read(fd, c, 1) == 1 ) && (*c != '\n') ) ;
      read(fd, c, 1);
    }
}

static void read_header_line(int fd, char *buf) { int i;
  i = 1;
  while( (read(fd, &buf[i], 1) == 1 ) && (buf[i] != '\n') && (buf[i] != '\r') && (i<79) ) i++;
  buf[i] = 0;
}

static int get_pgm_header(int fd, char *magic, int *width, int *height) { char buf[80];
  if( !match_key(fd, magic) ) return FALSE;
  read(fd, buf, 1);
  skip_comment(fd, '#', buf);
  read_header_line(fd, buf);
  sscanf(buf, "%d %d", width, height);
  read(fd, buf, 1);
  skip_comment(fd, '#', buf);
  read_header_line(fd, buf);
  return TRUE;
}

static int open_read(char *filename) { int fd;
  if( (fd = open(filename, O_BINARY|O_RDONLY)) < 0 )
    fprintf(stderr, "can't reset file `%s': %s\n", filename, strerror(errno));
  return fd;
}

static int open_read_pixmap(char *filename, char *magic, int *width, int *height) { int fd;
  if( (fd = open_read(filename)) < 0) return fd;
  if( !get_pgm_header(fd, magic, width, height) ) {
      fprintf(stderr, "can't read header of %s\n", filename);
      return -1;
    }
  return fd;
}

static unsigned char *alloc_pixmap(long size) { unsigned char *data;
  if( (data = (unsigned char *)malloc(size)) == NULL ) {
    fprintf(stderr, "malloc error\n");
    return NULL;
    }
  return data;
}

static void load_data(int fd, unsigned char *data, long size) { char *buffer;
  int count;

  buffer = (char *)data;
  while( size > 0 ) {
    count = IO_LEN;
    if( count > size ) count = size;
    read(fd, buffer, count);
    buffer += count;
    size -= count;
  }
}

static unsigned char *load_pixmap(char *filename, int *width, int *height) { int fd;
  long size;
  unsigned char *data;

  if( (fd = open_read_pixmap(filename, MAGIC_PGM, width, height)) < 0) return NULL;
  size = (long)*width * *height;
  data = alloc_pixmap(size);
  if( data != NULL ) load_data(fd, data, size);
  close(fd);
  return data;
}

static void put_header_line(int fd, char *buf) {
  write(fd, buf, strlen(buf));
}

static void put_header_info(int fd, char *mark, char *filename) { }

static void put_pgm_header(int fd, char *magic, int width, int height, char *filename) { char buf[80];
  put_header_line(fd, magic);
  put_header_info(fd, "# ", filename);
  sprintf(buf, "%d %d\n255\n", width, height);
  put_header_line(fd, buf);
}

static int open_write(char *filename) { int fd;
  if( (fd = open(filename, O_TRUNC|O_CREAT|O_BINARY|O_RDWR,0640)) < 0 )
/*
  if( (fd = open(filename, O_BINARY|O_CREAT|O_TRUNC|O_RDWR, S_IREAD|S_IWRITE)) < 0 )
*/
    fprintf(stderr, "can't Rewrite file `%s': %s\n", filename, strerror(errno));
  return fd;
}

static void store_data(int fd, unsigned char *data, long size) { char *buffer;
  int count;

  buffer = (char *)data;
  while( size > 0 ) {
    count = IO_LEN;
    if( count > size ) count = size;
    write(fd, buffer, count);
    buffer += count;
    size -= count;
    }
}

static void store_pixmap(char *filename, unsigned char *data, int width, int height) { int fd;
  if( (fd = open_write(filename)) < 0 ) return;
  put_pgm_header(fd, MAGIC_PGM, width, height, filename);
  store_data(fd, data, (long)width*height);
  close(fd);
}
        /*
                Lecture dnas le fichier "filename" de format pgm d'une image en niveau de gris  de rows lignes et cols colonnes
      		Le nombre de ligne et de colonnes est lu dans le fichier "filename"
		 L'image est cree dynamiquement, les pixels remplis puis l'image est retournee
                Retourne NULL en cas d'echec
        */
unsigned char ** lectureimagepgm(char *name, int* rows, int* cols) { unsigned char ** im;
  unsigned char * p;
  int i;
  p = load_pixmap(name,cols,rows);

  if (p==NULL) return NULL;
  if ( (im=calloc(*rows,sizeof(*im))) == NULL) return NULL;
  for (*im=p,i=1; i< *rows; i++) im[i] = im[i-1]+ *cols;
  return im;
}
        /*
                Ecriture dnas le fichier "filename" de format pgm d'une image en niveau de gris  de rows lignes et cols colonnes
        */
void ecritureimagepgm(char *name, unsigned char **im, int rows, int cols) {
  store_pixmap(name,*im,cols,rows);
}

        /*
                Ecriture dnas le fichier "filename"  d'une image de reels double precision
                Le format du fichier est "RAW", ie les donnees sont ecrites en binaire ligne par ligne et colonnes par colonnes.
        */
void ecritureimagedoubleraw(char *filename, double **im, int rows, int cols) { int fd;
  if( (fd = open_write(filename)) < 0 ) return;
  write(fd,*im,rows*cols*sizeof(**im));
  close(fd);
}

	/*
		Lecture et creation d'uen image de reels double precision contenu dns le fichier "filename"
		Le format du fichier est "RAW", ie les donnees sont ecrites en binaire ligne par ligne et colonnes par colonnes.
		Le nombre de ligne et de colonnes doivent etre connu
		L'image est cree dynamiquement, les pixels remplis
		et l'image est retournee
		Retourne NULL en cas d'echec
	*/
double** lectureimagedoubleraw( char *filename, int rows, int cols) { int fd;
  double** im=NULL;
  if( (fd = open_read(filename)) < 0 ) return NULL;
  if ( (im=alloue_image_double(rows,cols))==NULL) return NULL;
  read(fd,*im,rows*cols*sizeof(**im));
  close(fd);
}

        /*
                Creation d'une image de d'entiers 8bits non signes de nl lignes et nl colonnes
                Rq : l'allocation des donnees est faite en une seule allocation
                        on peut donc utiliser im[i][j] pour acceder a l'element d'indice i,j
                        ou bien *((*im)+i*nc+j) ou (*im)[i*nc+j] comme un tableau 1D
        */

unsigned char ** alloue_image(int nl, int nc) { int i;
  unsigned char** p;
  if ( (p=(unsigned char**)calloc(nl,sizeof(*p)))==NULL) return NULL;
  if ( (*p=(unsigned char*)calloc(nl*nc,sizeof(**p)))==NULL) return NULL;
  for (i=1; i<nl; i++) p[i]=p[i-1]+nc;
  return p;
}

	/*
		Creation d'une image de reelle double precision de nl lignes et nl colonnes
		Rq : l'allocation des donnees est faite en une seule allocation
			on peut donc utiliser im[i][j] pour acceder a l'element d'indice i,j
			ou bien *((*im)+i*nc+j) ou (*im)[i*nc+j] comme un tableau 1D
	*/
double ** alloue_image_double(int nl, int nc) { int i;
  double** p;
  if ( (p=(double **)calloc(nl,sizeof(*p)))==NULL) return NULL;
  if ( (*p=(double*)calloc(nl*nc,sizeof(**p)))==NULL) return NULL;
  for (i=1; i<nl; i++) p[i]=p[i-1]+nc;
  return p;
}

	/*
		Libere la memoire associe a l'image im
	*/
void libere_image(unsigned char** im) {
  free(*im); free(im);
}

	/*
		Conversion d'une image de 8bits non signes en une image de double
		Utilise en particulier pour fft
	*/
double** imuchar2double(unsigned char **im, int nl, int nc) { int i;
  double** res;
  if ( (res=alloue_image_double(nl,nc))==NULL) return NULL;
  for (i=0; i<nl*nc; i++) (*res)[i]=(double) (*im)[i];
  return res;
}

	/*
		Conversion d'une image de reels double en une image 8bits non signes
	*/
unsigned char**imdouble2uchar(double** im,int nl, int nc) { int i;
   unsigned char** res;
   if ( (res=alloue_image(nl,nc))==NULL) return NULL;
   for (i=0; i<nl*nc; i++) (*res)[i]=(unsigned char) ((*im)[i]+0.5);
   return res;
}

        /*
                Conversion d'une image de reels double en une image 8bits signes
        */
char**imdouble2char(double** im,int nl, int nc) { int i;
   char** res;
   if ( (res= (char**) alloue_image(nl,nc))==NULL) return NULL;
   for (i=0; i<nl*nc; i++) (*res)[i]=(char) ((*im)[i]+0.5);
   return res;
}

	/*
		Copie dasn une nouvelle image de la partie de l'image comprise entre les indices (oi,oj) et (fi,fj) de l'image im
        */
unsigned char** crop(unsigned char **im,int oi, int oj, int fi, int fj) { int i,j,nl,nc;
  unsigned char ** res;
  nl=fi-oi; nc=fj-oj;
  if ( (res=alloue_image(nl,nc))==NULL) return NULL;
  for(i=0; i<nl; i++)
     for(j=0; j<nc; j++)
       res[i][j]=im[oi+i][oj+j];
  return res;
}

/*
 Application du filtre gaussien directement sur l'image réelle et imaginaire.
 */
void gaussianFilter(double** reel, double** imag,
                        int nl, int nc, double ecartType) {
    double coef = -2*pow(M_PI*ecartType,2);
    for (int i=0; i<nl; ++i) {
        for (int j=0; j<nc; ++j) {
            double gauss = exp(coef*(pow((((double)i-nl/2)/nl),2)+pow((((double)j-nc/2)/nc),2)));
            reel[i][j] *= gauss;
            imag[i][j] *= gauss;
        }
    }
}

/**Application du filtre gaussien linéaire.
 *
 * @param reel
 *      Image d'origine sur laquelle appliquer le filtre.
 * @param nl
 *      Nombre de lignes de l'image.
 * @param nc
 *      Nombre de colonnes de l'image.
 * @param nlMask
 *      Nombre de lignes du masque.
 * @param ncMask
 *      Nombre de colonnes du masque.
 * @param ecartType
 *      Ecart type souhaité pour le filtre gaussien.
 * @return
 *      Image allouée contenant l'image d'origine filtrée.
 */
double** linearGaussianFilter(double** reel, int nl, int nc,
                            int nlMask, int ncMask, double ecartType) {
    //On initialise un tableau avec les exponentielles précalculées
    int max = (nlMask>ncMask)?nlMask:ncMask;
    double* exp_precalc = malloc(sizeof(double)*(max+1));
    for (int i=0; i<=max; ++i) {
        exp_precalc[i] = exp((double)-i*i/(2.0*pow(ecartType,2.0)));
    }
    //On calcule le résultat
    double** res = alloue_image_double(nl, nc);
    for (int i=0; i<nl; ++i) {
        for (int j=0; j<nc; ++j) {
            //On applique le masque sur ce pixel
            res[i][j] = 0.0;
            double normalize = 0.0;
            for (int u=-nlMask; u<=nlMask; ++u) {
                double somme = 0.0;
                for (int v=-ncMask; v<=ncMask; ++v) {
                    somme += exp_precalc[abs(v)]*reel[(i+u+nl)%nl][(j+v+nc)%nc];
                    normalize += exp_precalc[abs(v)]*exp_precalc[abs(u)];
                }
                res[i][j] += somme*exp_precalc[abs(u)];
            }
            res[i][j] /= normalize;
        }
    }
    free(exp_precalc);
    return res;
}

/**Application du filtre gaussien linéaire séparable.
 * L'image est modifiée directement.
 * @param reel
 *      Image d'origine sur laquelle appliquer le filtre.
 * @param nl
 *      Nombre de lignes de l'image.
 * @param nc
 *      Nombre de colonnes de l'image.
 * @param nlMask
 *      Nombre de lignes du masque.
 * @param ncMask
 *      Nombre de colonnes du masque.
 * @param ecartType
 *      Ecart type souhaité pour le filtre gaussien.
 */
void linearGaussianFilterSeparable(double** reel, int nl, int nc,
                            int nlMask, int ncMask, double ecartType) {
    //On initialise un tableau avec les exponentielles précalculées
    int max = (nlMask>ncMask)?nlMask:ncMask;
    double* exp_precalc = malloc(sizeof(double)*(max+1));
    for (int i=0; i<=max; ++i) {
        exp_precalc[i] = exp((double)-i*i/(2.0*pow(ecartType,2.0)));
    }
    //On calcule le résultat
    double** res = alloue_image_double(nl, nc);
    double normalize = 0.0;
    for (int u=-nlMask; u<=nlMask; ++u)
        normalize += exp_precalc[abs(u)];
    for (int i=0; i<nl; ++i) {
        for (int j=0; j<nc; ++j) {
            //On applique le masque sur ce pixel
            res[i][j] = 0.0;
            for (int u=-nlMask; u<=nlMask; ++u) {
                res[i][j] += reel[(i+u+nl)%nl][(j+nc)%nc]*exp_precalc[abs(u)];
            }
            res[i][j] /= normalize;
        }
    }
    normalize = 0.0;
    for (int v=-ncMask; v<=ncMask; ++v)
        normalize += exp_precalc[abs(v)];
    for (int i=0; i<nl; ++i) {
        for (int j=0; j<nc; ++j) {
            //On applique le masque sur ce pixel
            reel[i][j] = 0.0;
            for (int v=-ncMask; v<=ncMask; ++v) {
                reel[i][j] += exp_precalc[abs(v)]*res[(i+nl)%nl][(j+v+nc)%nc];
            }
            reel[i][j] /= normalize;
        }
    }
    free(*res);
    free(res);
    free(exp_precalc);
}

//----------------------Detection de contours-----------------------------------

/**Fonction permettant le débruitage d'une image
 * (à utilisée avant une detection de contours).
 *
 * @param entree
 *      Image à débruiter.
 * @param nlpoint
 *      Nombre de lignes de l'image.
 * @param ncpoint
 *      Nombre de colonnes de l'images.
 * @param filtrage
 *      0 pour appliquer un filtre gaussien, 1 pour appliquer un filtre laplacien.
 * @return
 *      Renvoie une image allouée contenant l'image d'origine filtrée.
 */
double** debruitageAux(  unsigned char** entree, int *nlpoint, int *ncpoint) {
  //L'argument filtrage permet de selectionner le type de filtre à appliquer:
  //0=filtre gaussien
  //1=filtre laplcien
  int nl=*nlpoint;
  int nc=*ncpoint;
  double **filtre=alloue_image_double(nl,nc);
  double pi = 3.141592654;
  double sigmax=3.0 ; //coefficient à ajuster en fonction de l'image source
  //convertion en double entree et filtre
  double **image=imuchar2double(entree,nl,nc);
  //padding avec des 0
  double **imagepadd=imuchar2double(entree,nl,nc);
  imagepadd=padimdforfft(image,&nl,&nc);
  //allocation de l'espace memoire
  double **fftImageRe=alloue_image_double(nl,nc);
  double **fftImageIm=alloue_image_double(nl,nc);
  double **zero=alloue_image_double(nl,nc);
  //fft image
  fft(imagepadd,zero,fftImageRe,fftImageIm,nl,nc);
  //shift image
  double **fftImageReshift=alloue_image_double(nl,nc);
  double **fftImageImshift=alloue_image_double(nl,nc);
  fftshift(fftImageRe,fftImageIm,fftImageReshift,fftImageImshift,nl,nc);
  //choix du filtre pour le filtrage
  gaussianFilter(fftImageReshift,fftImageImshift,nl,nc, sigmax);
  //shift
  double **shiftRe=alloue_image_double(nl,nc);
  double **shiftIm=alloue_image_double(nl,nc);
  fftshift(fftImageReshift,fftImageImshift,shiftRe,shiftIm,nl,nc);
  //fft inverse
  double **resultRe=alloue_image_double(nl,nc);
  double **resultIm=alloue_image_double(nl,nc);
  ifft(shiftRe,shiftIm,resultRe,resultIm,nl,nc);
  return resultRe;
}

/**Détection de contours du premier ordre (gradient)
 *
 * @param entree
 *      Image à débruiter.
 * @param nlpoint
 *      Nombre de lignes de l'image.
 * @param ncpoint
 *      Nombre de colonnes de l'images.
 * @param filtrage
 *      0 pour appliquer un filtre gaussien, 1 pour appliquer un filtre laplacien.
 * @return
 *      Renvoie une image allouée contenant l'image d'origine filtrée.
 */
double** ContourDetection(unsigned char** entree, int* nlpoint, int* ncpoint){
  int nc=*ncpoint;
  int nl=*nlpoint;
  double** Gx=alloue_image_double(nl,nc);
  double** Gy=alloue_image_double(nl,nc);
  double** module=alloue_image_double(nl,nc);
  double** filtre=alloue_image_double(3,3);
  double seuil=0; //Seuil de coupure pour le seuillage final
  //filtrage gaussien préalable
  double **image=debruitageAux(entree, nlpoint, ncpoint);
  //génération du filtre (ici Prewitt)
  for (int i=0; i<3;i++){
      filtre[i][0]=-1;
      filtre[i][1]=0;
      filtre[i][2]=1;
  }
  //calcul de la convolution de chaque pixel avec le filtre
  for (int x=1;x<nl-1;x++){ //On est pas sur les bords (nl-1)
    for (int y=1; y<nc-1;y++){
      double sumx=0;
      double sumy=0;
      for (int i=-1;i<=1;i++){
        for (int j=-1;j<=1;j++){
          sumx+=image[x+i][y+j]*filtre[i+1][j+1];
          sumy+=image[x+i][y+j]*filtre[j+1][i+1];
        }
      }
      Gx[x][y]=sumx;
      Gy[x][y]=sumy;
      //calcul du module
      module[x][y]=pow(sumx,2)+pow(sumy,2);
      //calcul de la valeur maximale, servant de référence pour le seuil
      if (pow(sumx,2)+pow(sumy,2)>seuil){
        seuil=pow(sumx,2)+pow(sumy,2);
      }
    }
  }
  for (int x=1;x<nl-1;x++){
    for (int y=1; y<nc-1;y++){
      //On coupe nos pixels à partir d'un certains pourcentage du seuil.
      //La valeur du pourcentage est à modifier en fonction de l'image source.
      if (module[x][y]>=(3.5/100.0)*seuil){
        module[x][y]=0; //Noir
      }else{
        module[x][y]=255;//Blanc
      }
    }
  }
  return module;
}

/**Détection de contours du second ordre (laplacien)
 *
 * @param entree
 *      Image à débruiter.
 * @param nlpoint
 *      Nombre de lignes de l'image.
 * @param ncpoint
 *      Nombre de colonnes de l'images.
 * @param filtrage
 *      0 pour appliquer un filtre gaussien, 1 pour appliquer un filtre laplacien.
 * @return
 *      Renvoie une image allouée contenant l'image d'origine filtrée.
 */
double** ContourDetectionLapl(unsigned char** entree, int* nlpoint, int* ncpoint){
  int nc=*ncpoint;
  int nl=*nlpoint;
  double** Lx=alloue_image_double(nl,nc);
  double** Ly=alloue_image_double(nl,nc);
  double** result=alloue_image_double(nl,nc);
  double** filtre=alloue_image_double(3,3);
  //filtrage gaussien préalable en plus du lissage du filtre
  double **image=debruitageAux(entree, nlpoint, ncpoint);
  //génération du filtre de laplace
  for (int i=0; i<3;i++){
    filtre[i][0]=1;
    filtre[i][1]=1;
    filtre[i][2]=1;
  }
  filtre[1][1]=-8;
  //Calcul de la convolution de chaque pixel avec le filtre
  for (int x=1;x<nl-1;x++){ //On est pas sur les bords
    for (int y=1; y<nc-1;y++){
      double sumx=0;
      double sumy=0;
      for (int i=-1;i<=1;i++){
        for (int j=-1;j<=1;j++){
          sumx+=image[x+i][y+j]*filtre[i+1][j+1];
          sumy+=image[x+i][y+j]*filtre[j+1][i+1];
        }
      }
      Lx[x][y]=sumx;
      Ly[x][y]=sumy;
    }
  }
  for (int x=1;x<nl-1;x++){
    for (int y=1; y<nc-1;y++){
      //Detection du saut du laplacien. On ne teste pas le signe mais la valeur par rapport
      //à 1 afin d'éliminer les bruits parasites
      if (Lx[x][y]*Lx[x+1][y] <=-1.0 || Ly[x][y]*Ly[x][y+1] <=-1.0){
        result[x][y]=0;
      }else{
        result[x][y]=255;
      }
    }
  }
  return result;
}
//-----------------------------------------------------------------------------------------

/**Enregistre la valeur du médian et enregistre la valeur dans l'argument median.
 *
 * @param histo
 *      Histogramme à partir duquel on calcule le médian
 * @param median
 *      Pointe vers l'ancienne valeur du médian, et doit être mis a jour
 * @param nb_left
 *      Pointe vers le nombre d'éléments situés à gauche (ou égal) du médian
 * @param moitie
 *      Le nombre minimum d'éléments à gauche à atteindre pour trouver le médian
 */
void findMedian(int8_t* histo, int* median, int* nb_left, int moitie) {
    while (*nb_left > moitie) {
        *nb_left -= histo[*median];
        (*median)--;
    }
    while (*nb_left < moitie) {
        (*median)++;
        *nb_left += histo[*median];
    }
}

/**Applique le filtre médian sur l'image "orig" en entrée, avec le
 * paramètre n (problèmes au bord pas pris en compte).
 *
 * @param orig
 *      Image d'origine sur laquelle appliquer le filtre.
 * @param nl
 *      Nombre de lignes de l'image.
 * @param nc
 *      Nombre de colonnes de l'image.
 * @param n
 *      Taille de la fenêtre.
 * @return
 *      Renvoie une image allouée contenant l'image d'origine filtrée.
 */
unsigned char** medianFilter(unsigned char** orig, int nl, int nc, int n) {
    unsigned char** res = alloue_image(nl,nc);
    int8_t histo[255] = {0};
    int old_median=0;
    int nb_left=0;
    int moitie = (2*n+1)*n+1;
    int i=n, j=n;
    int dir = 1;
    //Initialisation de l'histogramme
    for (int u=i-n; u<=i+n; u++) {
        for (int v=j-n; v<=j+n; v++) {
            histo[orig[u][v]]++;
        }
    }
    nb_left += histo[0];
    findMedian(histo, &old_median, &nb_left, moitie);
    res[i][j] = old_median;
    //On parcourt l'image en mettant à jour l'histogramme
    while (i < nc-n) {
        j += dir;
        if (j < n || j >= nl-n) {
            //Cette condition correspond à un changement de ligne  sur le parcourt de l'image.
            if (j >= nl-n-1) {
                j = nl-n-1;
                dir = -1;
                i++;
            }
            else {
                j = n;
                dir = 1;
                i++;
            }
            if (i >= nc-n) break;
            for (int k=j-n; k<=j+n; k++) {
                //On supprime les anciennes valeurs à l'histogramme
                int val_suppr = orig[i-n-1][k];
                histo[val_suppr]--;
                if (val_suppr <= old_median)
                    nb_left--;
                //On ajoute les nouvelles valeurs à l'histogramme
                int new_val = orig[i+n][k];
                histo[new_val]++;
                if (new_val <= old_median)
                    nb_left++;
            }
            //On cherche la valeur médiane en partant de l'ancienne
            findMedian(histo, &old_median, &nb_left, moitie);
            res[i][j] = old_median;
        }
        else {
            //On parcourt l'image sur toute sa ligne.
            for (int k=i-n; k<=i+n; k++) {
                //On supprime les anciennes valeurs à l'histogramme
                int val_suppr = orig[k][j-dir*(n+1)];
                histo[val_suppr]--;
                if (val_suppr <= old_median)
                    nb_left--;
                //On ajoute les nouvelles valeurs à l'histogramme
                int new_val = orig[k][j+dir*n];
                histo[new_val]++;
                if (new_val <= old_median)
                    nb_left++;
            }
            //On cherche la valeur médiane en partant de l'ancienne
            findMedian(histo, &old_median, &nb_left, moitie);
            res[i][j] = old_median;
        }
    }
    return res;
}

/**Applique le filtre bilatéral sur l'image "orig" en entrée, avec les
 * paramètres sigma1 et sigma2 (problèmes au bord pas pris en compte).
 *
 * @param orig
 *      Image d'origine sur laquelle appliquer le filtre.
 * @param nl
 *      Nombre de lignes de l'image.
 * @param nc
 *      Nombre de colonnes de l'image.
 * @param sigma1
 *      Paramètre lié à la taille des détails à préserver.
 * @param sigma2
 *      Paramètre lié au bruit de l'image.
 * @return
 *      Renvoie une image allouée contenant l'image d'origine filtrée.
 */
double** bilateralFilter(double** orig,
                        int nl, int nc, int sigma1, double sigma2) {
    double** res = alloue_image_double(nl,nc);
    for (int x=3*sigma1; x<nl-3*sigma1; ++x) {
        for (int y=3*sigma1; y<nc-3*sigma1; ++y) {
            double num = 0.0;
            double denom = 0.0;
            for (int i=-3*sigma1; i<=3*sigma1; ++i) {
                for (int j=-3*sigma1; j<=3*sigma1; ++j) {
                    double coef = exp(-(double)(i*i+j*j)/(2.0*sigma1*sigma1))*
                            exp(-(double)pow(orig[x+i][y+j]-orig[x][y],2.0)/(2.0*sigma2*sigma2));
                    num += coef*orig[x+i][y+j];
                    denom += coef;
                }
            }
            res[x][y] = num/denom;
        }
    }
    return res;
}
