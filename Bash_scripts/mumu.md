# Mumu installation and usage on the cluster

## installation

It basically follow the [README](https://github.com/frederic-mahe/mumu) of mumu

Get the .tar.gz of the desired version and untar it (here I do it in the ~/work/app/ folder):

```bash
wget wget https://github.com/frederic-mahe/mumu/archive/refs/tags/v1.1.4.tar.gz -O mumu_v1.1.4.tar.gz

tar -xvf mumu_v1.1.4.tar.gz
```

Now go to the mumu folder :

```bash
cd mumu-1.1.4/
```

Load the gcc compiler :

```bash
module load compilers/gcc/15.1.0
```

Then follow the steps of the readme :

```bash
make
make check
```

But then the install cannot be done using make install because it needs root privilege.
So we just export the paths to the the bash environnement :

```bash
export PATH="$HOME/work/app/mumu-1.1.4:$PATH"
```

Then we reload the environment

```bash
source ~/.bashrc
```

Test it :

```bash
mumu -h
```

bonus : adding mumu to the man

First create a man folder :

```bash
mkdir -p ~/.local/share/man/man1
```

Then add the man file (here doing the command from the mumu folder)

```bash
cp man/mumu.1 ~/.local/share/man/man1/
```

Add the new man folder to the MANPATH

```bash
export MANPATH="$HOME/.local/share/man:$MANPATH"
source ~/.bashrc
```

Then the final test

```bash
man mumu
```

It works ! (tested on 09/08/2026)
