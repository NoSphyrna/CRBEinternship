# Installation of the TaxInfo package on the Genotoul cluster

## Initialization

Once connected to the cluster, create a folder R to contain the packages:

```bash
# Create a directory named 'R' elsewhere (here in ~/work)
mkdir ~/work/R
# Create a symbolic link to this new directory
ln -sr ~/work/R ~/R
```

## Session launch and loading of modules (with slurm)

Then launch a session with enough ram to avoid an OOM (Out Of Memory)
during the installation of the packages

(Here we allocate 8Go of RAM for our installation session)

```bash
srun --pty --mem=8G bash
```

And if we wish for it to be a bit faster, we can allocate more cpus:

```bash
srun --pty --mem=8G --cpus-per-task=4 bash
```

Launch the modules of the lastest versions of R and gcc
(gcc to be able to compile the package "sf" which contains
"abseil-cpp-devel" which depends on TaxInfo)

```bash
module load compilers/gcc/15.1.0 
module load statistics/R/4.5.0
```

If needed, check the lastest versions installed on the cluster with:

```bash
module avail
```

Or more specifically the modules gcc and R :

```bash
module avail compilers/gcc 
module avail statistics/R
```

## Installation of abseil : (necessary for the package s2 itself necessary for taxinfo but not for MiscMetabar)

We create first the directory where we wish to store the package abseil in a
permanent way and the directory where we wish to do the compilation of the package.

Here for example, I used the directory **~/work/local** for the permanent storage
and **~/work/tmp** for the temporary storage for the compilation.

```bash
mkdir ~/work/local
mkdir ~/work/tmp
```

We place ourselves now in the temporary file:

```bash
cd ~/work/tmp
```

We get the most recent version of abseil from the
[abseil github archives](https://github.com/abseil/abseil-cpp/releases):

```bash
wget https://github.com/abseil/abseil-cpp/releases/download/20260107.1/abseil-cpp-20260107.1.tar.gz
tar xvf abseil-cpp-20260107.1.tar.gz
```

We are in the abseil-cpp directory which has just been created
and we create the build directory and we enter it:

```bash
cd abseil-cpp-20260107.1
mkdir build && cd build
```

Now, it is necessary to configure the compilation with cmake:

Explanation of the command:

- **cmake .. \\**: we launch cmake from the parent folder
  (".." = "abseil-cpp-20260107.1" here), the "\\" allows to write the command
  on several lines for more readability but we can also write everything on
  a single line (see below).
- **-DCMAKE_INSTALL_PREFIX=~/work/local**: The folder where we wish to
  install the compiler files.
- **-DCMAKE_CXX_STANDARD=17**: The version of the gcc compiler minimal for abseil
  (allows to check at the compilation that we have the right version of gcc)
- **-DABSL_PROPAGATE_CXX_STD=ON**: Propagation of the previous verification to
  all internal modules.
- **-DCMAKE_BUILD_TYPE=Release**: Compression and optimization of the binaries
  created by the compilation so that the package is lighter.
- **-DBUILD_SHARED_LIBS=ON**: Creation of dynamic libraries (.so) instead of
  static ones so that R can load them.

```bash
cmake .. \
  -DCMAKE_INSTALL_PREFIX=~/work/local \
  -DCMAKE_CXX_STANDARD=17 \
  -DABSL_PROPAGATE_CXX_STD=ON \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON
```

The version in one line of this same command:

```bash
cmake .. -DCMAKE_INSTALL_PREFIX=~/work/local -DCMAKE_CXX_STANDARD=17 -DABSL_PROPAGATE_CXX_STD=ON -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON
```

Now, we need to compile the package,
**it is necessary to choose ONE command from the 2 following :**

If we have chosen to add CPUs at the
[launch of the session](<#Session_launch_and_loading_of_modules_(with_slurm)>)
to go faster we can execute the following command with the same number
of CPUs that we have allocated.

```bash
make -j4
```

Otherwise, we simply execute:

```bash
make
```

Finally, we install the compiled files in the folder **~/work/local**
previously created :

```bash
make install
```

We can clean the temporary folder:

```bash
cd ~/work/tmp
rm -rf abseil-cpp-20260107.1 abseil-cpp-20260107.1.tar.gz
```

Finally, it is necessary to export the environment variables so that the R
libraries know where to find abseil. And so that it is permanent, we write
the export commands in the file ".bashrc". Then we reload the file
".bashrc" with "source".

```bash
echo 'export PKG_CONFIG_PATH=$HOME/work/local/lib64/pkgconfig:$PKG_CONFIG_PATH' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=$HOME/work/local/lib64:$LD_LIBRARY_PATH' >> ~/.bashrc
source ~/.bashrc 
```

## Installation of the packages in R

Launch R :

```bash
R
```

Then install the package "sf" :

```R
install.packages("sf")
```

Install package "s2" :

```R
install.packages("s2")
```

Install MiscMetabar :

```R
pak::pak("adrientaudiere/MiscMetabar")
```

And finally install TaxInfo :

```R
pak::pak("adrientaudiere/taxinfo")
```

> Note: Here the package should install itself despite the package "abseil-cpp-devel"
> absent. Indeed, we have installed manually abseil but pak does not notice it.
> So this warning is normal.

## Example of a complete functional installation (tried the 20/05/2026)

```bash
# Creation of the R directory for
mkdir ~/work/R
# Create a symbolic link to this new directory
ln -sr ~/work/R ~/R

# Launch of a session and loading of the modules

srun --pty --mem=8G --cpus-per-task=4 bash
module load compilers/gcc/15.1.0 
module load statistics/R/4.5.0

# Installation of abseil

## Creation of the compilation files and download
mkdir ~/work/local
mkdir ~/work/tmp

cd ~/work/tmp

##  Retrieval of the lastest version of abseil
wget https://github.com/abseil/abseil-cpp/releases/download/20260107.1/abseil-cpp-20260107.1.tar.gz
tar xvf abseil-cpp-20260107.1.tar.gz

## Creation of the compilation file
cd abseil-cpp-20260107.1
mkdir build && cd build

## Configuration of the compilation
cmake .. \
  -DCMAKE_INSTALL_PREFIX=~/work/local \
  -DCMAKE_CXX_STANDARD=17 \
  -DABSL_PROPAGATE_CXX_STD=ON \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_SHARED_LIBS=ON

## Compilation and exportation of the packages towards ~/work/local/
make -j4
make install

# Cleaning of the tmp folder
cd ~/work/tmp
rm -rf abseil-cpp-20260107.1 abseil-cpp-20260107.1.tar.gz

# We return to the home
cd 

# Export of the environment variables
echo 'export PKG_CONFIG_PATH=$HOME/work/local/lib64/pkgconfig:$PKG_CONFIG_PATH' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=$HOME/work/local/lib64:$LD_LIBRARY_PATH' >> ~/.bashrc
source ~/.bashrc 

# Launch of R
R
```

```R
# Installation of packages
install.packages("sf")

install.packages("s2")

pak::pak("adrientaudiere/MiscMetabar")

pak::pak("adrientaudiere/taxinfo")
```
