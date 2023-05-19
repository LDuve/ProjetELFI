femFullSystemConstrain

Cette fonction permet d'appliquer une contrainte de Dirichlet (condition de bord) sur un système linéaire. Elle fixe une valeur spécifiée pour un nœud donné dans le système,
 en modifiant les matrices A et B. Les autres éléments de la ligne et de la colonne correspondant au nœud sont mis à zéro.

Arguments :
mySystem : Pointeur vers la structure femFullSystem représentant le système linéaire.
myNode : Index du nœud auquel la contrainte sera appliquée.
myValue : La valeur à laquelle le nœud spécifié sera contraint.




femFullSystemConstrainNeumann
Cette fonction permet d'appliquer une condition de Neumann (flux) sur un système linéaire. Elle ajoute une valeur spécifiée à l'élément correspondant au nœud dans le vecteur B.

Arguments :
mySystem : Pointeur vers la structure femFullSystem représentant le système linéaire.
myNode : Index du nœud auquel la condition de Neumann sera appliquée.
myValue : La valeur à ajouter au nœud spécifié dans la condition de Neumann.




femElasticityAddBoundaryCondition
Cette fonction permet d'ajouter une condition de bord à un problème d'élasticité. 
Elle prend en compte le nom du domaine, le type de condition de bord (Dirichlet ou Neumann) et la valeur associée. 
Elle met à jour la liste des conditions de bord du problème, ainsi que les nœuds contraints.

Arguments :
theProblem : Pointeur vers la structure femProblem représentant le problème d'élasticité.
nameDomain : Nom du domaine auquel la condition aux limites sera appliquée.
type : Type de la condition aux limites (DIRICHLET_X, DIRICHLET_Y, NEUMANN_X, NEUMANN_Y).
value : La valeur associée à la condition aux limites.




femElasticitySolveSym
Cette fonction résout un problème d'élasticité en utilisant une approche symétrique.

Arguments: 
theProblem : Pointeur vers la structure femProblem représentant le problème d'élasticité.
solver : Indicateur du solveur à utiliser (1 pour Gauss, 2 pour Bande, 3 pour Frontal).
Sortie : 
La fonction renvoie un pointeur vers un tableau de valeurs doubles représentant la solution du système.




femFullSystemEliminateBande
Cette fonction résout un système linéaire à l'aide de la méthode d'élimination bande.

Arguments :
mySystem : Pointeur vers la structure femFullSystem représentant le système linéaire.
Sortie : 
La fonction renvoie un pointeur vers un tableau de valeurs doubles représentant la solution du système.




femFullSystemEliminate
Cette fonction résout un système linéaire à l'aide de la méthode d'élimination gaussienne.

Arguments :
mySystem : Pointeur vers la structure femFullSystem représentant le système linéaire.
Sortie :
La fonction renvoie un pointeur vers un tableau de valeurs doubles représentant la solution du système.




femFullSystemEliminateFrontal
Cette fonction résout un système linéaire à l'aide de la méthode d'élimination frontale.

Arguments :
mySystem : Pointeur vers la structure femFullSystem représentant le système linéaire.
Sortie :
La fonction renvoie un pointeur vers un tableau de valeurs doubles représentant la solution du système.




calculateBandwidth
Cette fonction calcule la largeur de la bande d'une matrice donnée.

Arguments : 
A : Matrice représentée sous forme d'un tableau 2D de valeurs doubles.
size : Taille de la matrice (nombre de lignes et de colonnes).
Sortie : 
La fonction renvoie un entier représentant la largeur de la bande de la matrice.