# CRBEinternship

Main scripts for the internship at CRBE

## Installation du package TaxInfo sur le cluster Genotoul

### Initialisation

Une fois connecté sur le cluster, créer un dossier R pour
contenir les packages :

```bash
# Create a directory named 'R' elsewhere (here in ~/work)
mkdir ~/work/R
# Create a symbolic link to this new directory
ln -sr ~/work/R ~/R
```

### Lancement de session et chargement des modules

Ensuite lancer une session avec suffisament de ram pour évitier
un OOM (Out Of Memory) lors de l'installation des packages

(Ici on alloue 8Go de RAM pour notre session d'installation)

```bash
srun --pty --mem=8G bash
```

Et si l'on souhaite que ce soit un peu plus rapide, on peut allouer plus de cpus:

```bash
srun --pty --mem=8G --cpus-per-task=4 bash
```

Lancer les modules des dernières version de R et de gcc (gcc
pour pouvoir compiler le package "sf" qui contient "abseil-cpp-devel"
dont dépend TaxInfo)

```bash
module load compilers/gcc/15.1.0 
module load statistics/R/4.5.0
```

Si besoin, vérifier les dernières versions installées sur le
cluster avec :

```bash
module avail
```

Ou plus spécifiquement les modules gcc et R :

```bash
module avail compilers/gcc 
module avail statistics/R
```

### Installation d'abseil : (nécessaire pour le package s2 lui-même nécessaire pour taxinfo mais pas pour MiscMetabar)

On crée d'abord l'endroit où l'on souhaite stocker le package abseil de
manière permanante et l'endroit ou l'on souhaite faire la compilation du package.

Ici par example, j'ai utiliser le répertoire **~/work/local** pour le stockage permanant
et **~/work/tmp** pour le stockage temporaire pour la compilation.

```bash
mkdir ~/work/local
mkdir ~/work/tmp
```

On se place maintenant dans le fichier temporaire :

```bash
cd ~/work/tmp
```

On récupère la version la plus récente d'abseil à partir des
[archives github d'abseil](https://github.com/abseil/abseil-cpp/releases) :

```bash
wget https://github.com/abseil/abseil-cpp/releases/download/20260107.1/abseil-cpp-20260107.1.tar.gz
tar xvf abseil-cpp-20260107.1.tar.gz
```

On dans le répertoire abseil-cpp qui bient d'être créé
et on crée le répertoire de build et on y rentre :

```bash
cd abseil-cpp-20260107.1
mkdir build && cd build
```

Maintenant, il faut configurer la compilation avec cmake :

Explication de la commande :

- **cmake .. \\** : on lance cmake depuis le dossier parent
  (".." = "abseil-cpp-20260107.1" ici), le "\\" permet d'écrire la commande
  sur plusieurs lignes pour plus de lisibilité mais on peut aussi
  tout écrire sur une seule ligne (voir en dessous).
- **-DCMAKE_INSTALL_PREFIX=~/work/local** : Le dossier où l'on souhaite
  installer les fichier compiler.
- **-DCMAKE_CXX_STANDARD=17** : La version du compilateur gcc minimale pour
  abseil (permet de vérifier à la compilation que l'on a bien la bonne version
  de gcc)
- **-DABSL_PROPAGATE_CXX_STD=ON** : Propagation de la vérification précédente
  à tous les modules internes.
- **-DCMAKE_BUILD_TYPE=Release** : Compression et optimisation des binaires créés
  par la compilation pour que le package soit plus léger.
- **-DBUILD_SHARED_LIBS=ON** : Création de bibliotèques dynamiques (.so) au lieu
  de statiques pour que R puisse les charger.

```bash
cmake .. \
  -DCMAKE_INSTALL_PREFIX=~/work/local \
  -DCMAKE_CXX_STANDARD=17 \
  -DABSL_PROPAGATE_CXX_STD=ON \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON
```

La version en une ligne de cette même commande:

```bash
cmake .. -DCMAKE_INSTALL_PREFIX=~/work/local -DCMAKE_CXX_STANDARD=17 -DABSL_PROPAGATE_CXX_STD=ON -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON
```

Maintenant, il faut compiler le package, **il faut choisir UNE commande parmis les 2 suivantes :**

Si on a choisi d'ajouter des CPUs à lors du
[lancement de la session](#lancement-de-session-et-chargement-des-modules)
pour aller plus vite on peut exécuter la commande suivante avec le même nombre
de CPUs que l'on a alloué.

```bash
make -j4
```

Sinon, il faut simplement exécuter :

```bash
make
```

Enfin, on installe les fichier compilés dans le dossier **~/work/local**
précédemment créé :

```bash
make install
```

On peut nettoyer le dossier temporaire :

```bash
cd ~/work/tmp
rm -rf abseil-cpp-20260107.1 abseil-cpp-20260107.1.tar.gz
```

Enfin, il faut exporter les variables d'environnement pour que les librairies R
sachent où trouver abseil.
Et pour que ça soit permanant, on écrit les commandes d'export dans le fichier ".bashrc".
Puis on recharge le fichier ".bashrc" avec "source".

```bash
echo 'export PKG_CONFIG_PATH=$HOME/work/local/lib64/pkgconfig:$PKG_CONFIG_PATH' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=$HOME/work/local/lib64:$LD_LIBRARY_PATH' >> ~/.bashrc
source ~/.bashrc 
```

### Installation des packages dans R

Lancer R :

```bash
R
```

Puis installer le package "sf" :

```R
install.packages("sf")
```

Installer le package "s2" :

```R
install.packages("s2")
```

Installer MiscMetabar :

```R
pak::pak("adrientaudiere/MiscMetabar")
```

Et infin installer TaxInfo :

```R
pak::pak("adrientaudiere/taxinfo")
```

> Note: Ici le package devrait s'installer malgré le package "abseil-cpp-devel"
> absent. En effet, on a installé manuellement abseil mais pak ne le remarque
> pas. Donc ce warning est normal.

### Exemple d'installation complète fonctionnelle au 20/05/2026

```bash
# Création du répertoir R pour 
mkdir ~/work/R
# Create a symbolic link to this new directory
ln -sr ~/work/R ~/R

# Lancement d'un session et chargement des modules

srun --pty --mem=8G --cpus-per-task=4 bash
module load compilers/gcc/15.1.0 
module load statistics/R/4.5.0

# Installation d'abseil

## Creation des fichiers de compilation et téléchargement
mkdir ~/work/local
mkdir ~/work/tmp

cd ~/work/tmp

## Récupération de la dernière version d'abseil
wget https://github.com/abseil/abseil-cpp/releases/download/20260107.1/abseil-cpp-20260107.1.tar.gz
tar xvf abseil-cpp-20260107.1.tar.gz

## création du fichier de compilation
cd abseil-cpp-20260107.1
mkdir build && cd build

## Configuration de la compilation
cmake .. \
  -DCMAKE_INSTALL_PREFIX=~/work/local \
  -DCMAKE_CXX_STANDARD=17 \
  -DABSL_PROPAGATE_CXX_STD=ON \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON

## Compilation et exportation des paquets vers ~/work/local/
make -j4
make install

# Nettoyage du dossier tmp
cd ~/work/tmp
rm -rf abseil-cpp-20260107.1 abseil-cpp-20260107.1.tar.gz

# On revient au home
cd 

# Export des variables d'envirronement
echo 'export PKG_CONFIG_PATH=$HOME/work/local/lib64/pkgconfig:$PKG_CONFIG_PATH' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=$HOME/work/local/lib64:$LD_LIBRARY_PATH' >> ~/.bashrc
source ~/.bashrc 

# Lancement de R
R
```

```R
# Installation des packages
install.packages("sf")

install.packages("s2")

pak::pak("adrientaudiere/MiscMetabar")

pak::pak("adrientaudiere/taxinfo")
```
