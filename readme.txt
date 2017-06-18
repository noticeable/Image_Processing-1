MANUEL D'UTILISATION DU PROGRAMME

1)FILTRE GAUSSIEN ET FILTRE GAUSSIEN LINEAIRE

Le filtre gaussien et filtre gaussien linéaire ont été implémentés dans pgm.c.
Le programme c 'part1.c' utilise ces filtres pour comparer le PSNR fourni ou enregistrer l'image filtrée,
et donne ainsi un exemple d'utilisation de ces filtres.
Pour le compiler :
	make part1

Ensuite, il suffit de l'exécuter.

2)DETECTION DE CONTOURS

Afin d'utiliser les programmes de détection de contours, il faut dans un
premier temps compiler ces derniers comme suit:
    make contours
    make contourslaplacien

Une fois cette opération effectuée, il suffit de lancer ces exécutable en spécifiant
l'image source comme premier argument et l'image destination comme second argument.

Afin de pouvoir améliorer la détection, l'utilisateur peut modifier certains paramètres du programme:
    -Dans la fonction debruitageAux, le paramètre sigmax permet de regler l'intensitée du filtrage
    -Dans la fonction ContoursDetection, le poucentage de seuillage est également modifiable.
