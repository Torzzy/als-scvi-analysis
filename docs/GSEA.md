# Explication de GSEA

Pour calculer les pathways enrichis, GSEA prend la liste des gènes triée d'une manière quelconque (dans le sens où plusieurs rankings peuvent être pertinents).
Pour chaque pathway de sa base de données, elle va calculer un score d'enrichissement de la manière suivante :

Supposons que G<sub>pathway</sub> soit l'ensemble des gènes d'un pathway et G<sub>DE</sub> l'ensemble des gènes de notre étude.

On va calculer une somme partielle de la manière suivante :

## Calcul de l'ES (Enrichment Score) :

<img src="../sources/latex/ES.png" width="300">

On parcourt chaque gène de notre dataset. S'il est dans le pathway que l'on étudie, alors on ajoute la quantité de la première ligne. Cette quantité favorise les gènes bien classés si on pose p=2 par exemple.
Si le gène n'y est pas, alors on retire une petite quantité dépendant de la taille de G<sub>pathway</sub>.

On parcourt comme ça tous les gènes dans l'ordre trié et on ajoute le score du gène i dans la somme partielle à chaque itération.
Pendant ce calcul, on garde en mémoire la valeur minimale et maximale qu'a prise cette somme partielle. 

L'idée étant que les pathways ayant le plus petit ou le plus grand score parmi tous les pathways soient ceux qui nous intéressent.

## Calcul du NES (Normalized Enrichment Score) :

On définit le NES pour un pathway de la manière suivante :

<img src="../sources/latex/NES.png" width="200">

L'idée est de comparer les scores obtenus pour chaque pathway à ce qu'on pourrait obtenir par hasard.
Pour cela, on fait un certain nombre de tirages aléatoires de rankings (1000 dans notre cas) et on regarde en moyenne quel est le score du pathway.
Le NES correspond au ratio entre le score de notre pathway avec notre ranking et le score moyen avec un ranking au hasard.

Ce calcul permet d'obtenir une p-value empirique puis de calculer une FDR par pathway.