#include "pgm.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>

/** file est le répertoire dans lequel se trouve l'image (doit contenir un /)
 * L'image filtrée est alors enregistrée dans le sous dossier 'subDirectoty'
 * dans le premier dossier de file, en conservant la suite du chemin de file
 * (Ce sous-dossier doit exister).
 */
char * createPath(const char *file, const char *subDirectory, const double e) {
    // 5 en plus pour "/" et le nombre _e
    char * path = calloc(sizeof(char), strlen(file)+strlen(subDirectory)+6);
    int i=0;
    while (file[i] != '/' && file[i] != '\0') {
        path[i] = file[i];
        i++;
    }
    if (file[i] == '\0') {
        printf("Un / est attendu dans le chemin vers image bruitée");
        exit(1);
    }
    path[i++] = '/';
    int save = i;
    strcpy(&(path[i]), subDirectory);
    i += strlen(subDirectory);
    //Si le sous dossier n'existe pas on le crée.
    struct stat st = {0};
    if (stat(path, &st) == -1) {
        mkdir(path, 0700);
    }
    path[i++] = '/';
    strcpy(&(path[i]), &(file[save]));
    i += strlen(&(file[save]))-4;
    save = strlen(file)-4;
    //On ajoute le nombre e
    path[i++] = '_';
    char nb[4];
    sprintf(nb, "%3.0f", e*100.0);
    for (int j=0; j<3; j++)
        if (nb[j] == ' ')
            nb[j] = '0';
    strcpy(&(path[i]), nb);
    i += 3;
    strcat(&(path[i]), &(file[save]));
    return path;
}

/** file_orig représente le fichier image sans bruit, et
 * file_bruit représente le fichier image avec bruit.
 */
double comparePSNRGaussianFilter(char* file_orig, char* file_bruit, double ecart_type) {
    /****** Chargement des images *******/
    int nl1,nc1,nl2,nc2;
    unsigned char** img_orig = lectureimagepgm(file_orig,&nl1,&nc1);
    if (img_orig==NULL)  { puts("Lecture image impossible"); exit(1); }
    unsigned char** img_bruit = lectureimagepgm(file_bruit,&nl2,&nc2);
    if (img_bruit==NULL)  { puts("Lecture image impossible"); exit(1); }
    assert(nl1==nl2 && nc1==nc2);
    /******* Application du filtre Gaussien ********/
    double** im1 = imuchar2double(img_bruit,nl1,nc1);
    im1 = padimdforfft(im1,&nl1,&nc1);
    double** im2 = alloue_image_double(nl1,nc1);
    double** im3 = alloue_image_double(nl1,nc1);
    double** im4 = alloue_image_double(nl1,nc1);
    fft(im1, im2, im3, im4, nl1, nc1);
    fftshift(im3, im4, im1, im2, nl1, nc1);
    gaussianFilter(im1, im2, nl1, nc1, ecart_type);
    fftshift(im1, im2, im3, im4, nl1, nc1);
    ifft(im3, im4, im1, im2, nl1, nc1);
    unsigned char** img_filtre = crop(imdouble2uchar(im1,nl1,nc1),0,0,nl2,nc2);
    /******* On enregistre l'image filtree pour voir le résultat *******/
    //ecritureimagepgm(createPath(file_bruit, "gaussianFilter", ecart_type), img_filtre, nl2, nc2);
    /********* Calcule de l'efficacité du filtre avec PSNR **********/
    return psnr(img_orig, img_filtre, nl2, nc2);
}

double comparePSNRLinearGaussianFilter(char* file_orig, char* file_bruit,
                                       int nlMask, int ncMask, double ecart_type) {
    int nlOrig, ncOrig;
    unsigned char** img_orig = lectureimagepgm(file_orig,&nlOrig,&ncOrig);
    if (img_orig==NULL)  { puts("Lecture image impossible"); exit(1); }
    int nl, nc;
    unsigned char** im1=lectureimagepgm(file_bruit,&nl,&nc);
    if (im1==NULL)  { puts("Lecture image impossible"); exit(1); }
    assert(nlOrig==nl && ncOrig==nc);
    /******* Application du filtre Gaussien linéaire********/
    double** im2 = imuchar2double(im1,nl,nc);
    linearGaussianFilterSeparable(im2, nl, nc, nlMask, ncMask, ecart_type);
    unsigned char** img_filtre = imdouble2char(im2, nl, nc);
    /******* On enregistre l'image filtree pour voir le résultat *******/
    //ecritureimagepgm(createPath(file_bruit, "linearGaussianFilter", ecart_type), img_filtre, nl, nc);
    /********* Calcule de l'efficacité du filtre avec PSNR **********/
    return psnr(img_orig, img_filtre, nl, nc);
}

int main (int ac, char **av) {
    /** Compare différentes valeurs de PSNR **/
    double res;
    double ecartType = 1.0;
    printf("Comparaison formes2 et formes2bb10 :\n");
    res = comparePSNRGaussianFilter("imgs/formes2.pgm", "imgs/formes2bb10.pgm", ecartType);
    printf("PSNR = %f\n", res);
    res = comparePSNRLinearGaussianFilter("imgs/formes2.pgm", "imgs/formes2bb10.pgm",
                                          4*ecartType,4*ecartType,ecartType);
    printf("PSNR (linear) = %f\n", res);

    printf("Comparaison formes2 et formes2pets1 :\n");
    res = comparePSNRGaussianFilter("imgs/formes2.pgm", "imgs/formes2pets1.pgm", ecartType);
    printf("PSNR = %f\n", res);
    res = comparePSNRLinearGaussianFilter("imgs/formes2.pgm", "imgs/formes2pets1.pgm",
                                          4*ecartType,4*ecartType,ecartType);
    printf("PSNR (linear) = %f\n", res);

    printf("Comparaison formes2 et formes2sp1 :\n");
    res = comparePSNRGaussianFilter("imgs/formes2.pgm", "imgs/formes2sp1.pgm", ecartType);
    printf("PSNR = %f\n", res);
    res = comparePSNRLinearGaussianFilter("imgs/formes2.pgm", "imgs/formes2sp1.pgm",
                                          4*ecartType,4*ecartType,ecartType);
    printf("PSNR (linear) = %f\n", res);

    /** Enregistre différentes valeurs de PSNR dans un fichier **/
//    FILE* fichier = NULL;
//    fichier = fopen("data.txt", "w");
//    if (fichier == NULL) exit(1);
//    fprintf(fichier, "Gaussien : ");
//    fprintf(fichier, "%f ", comparePSNRGaussianFilter("imgs/formes2.pgm", "imgs/formes2bb50.pgm", ecartType));
//    fprintf(fichier, "\nLinear Gaussien : ");
//    for (int i=1; i<=6*ecartType; i++) {
//        fprintf(fichier, "%f ",
//                comparePSNRLinearGaussianFilter("imgs/formes2.pgm", "imgs/formes2bb50.pgm", i, i, ecartType));
//    }
//    fclose(fichier);

    /** Calcule des différences en temps de calcule **/
//    FILE* fichier = NULL;
//    fichier = fopen("time.txt", "w");
//    if (fichier == NULL) exit(1);
//    double sigma = 0.0;
//    int nl = 256, nc = 256;
//    double** im1 = alloue_image_double(nl,nc);
//    double** im2 = alloue_image_double(nl,nc);
//    double** im3 = alloue_image_double(nl,nc);
//    double** im4 = alloue_image_double(nl,nc);
//    fprintf(fichier, "\nTime gauss :\n");
//    for (sigma = 0.2; sigma<6.0; sigma+=0.2) {
//        double gauss = 0.0;
//        clock_t debut, fin;
//        for (int i=0; i<10; i++) {
//            debut = clock();
//            fft(im1, im2, im3, im4, nl, nc);
//            fftshift(im3, im4, im1, im2, nl, nc);
//            gaussianFilter(im1, im2, nl, nc, sigma);
//            fftshift(im1, im2, im3, im4, nl, nc);
//            ifft(im3, im4, im1, im2, nl, nc);
//            fin = clock();
//            gauss += ((double)fin-debut)/CLOCKS_PER_SEC;
//        }
//        gauss /= 10.0;
//        fprintf(fichier, "%f ", gauss);
//    }
//    fprintf(fichier, "\nTime linear :\n");
//    for (sigma = 0.2; sigma<6.0; sigma+=0.2) {
//        double linear = 0.0;
//        clock_t debut, fin;
//        for (int i=0; i<10; i++) {
//            debut = clock();
//            linearGaussianFilter(im1, nl, nc, 3.5*sigma, 3.5*sigma, sigma);
//            fin = clock();
//            linear += ((double)fin-debut)/CLOCKS_PER_SEC;
//        }
//        linear /= 10.0;
//        fprintf(fichier, "%f ", linear);
//    }
//    fprintf(fichier, "\nTime separable :\n");
//    for (sigma = 0.2; sigma<6.0; sigma+=0.2) {
//        double linear = 0.0;
//        clock_t debut, fin;
//        for (int i=0; i<10; i++) {
//            debut = clock();
//            linearGaussianFilterSeparable(im1, nl, nc, 3.5*sigma, 3.5*sigma, sigma);
//            fin = clock();
//            linear += ((double)fin-debut)/CLOCKS_PER_SEC;
//        }
//        linear /= 10.0;
//        fprintf(fichier, "%f ", linear);
//    }
//    fclose(fichier);
    return 0;
}
