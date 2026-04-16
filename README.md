# cg-dsDAN-openmm

This program is initially prepared to implement the double strand DNA model under the framework of HPS

Author info: Qiang Zhu (Vice-Chancellor Research Fellow at University of Wollongong)

## Architecture

![architecture drawn with drawio](./README/HPS_CG.drawio.png)


## How to install the environment

The program runs with `python3`. The philosophy is to minimize the dependence of external libraries. 

Two mainly external libraries, namely, `OpenMM` and `Biopython` are required.

For your convenience, we also provide you with a environment file [`requirements.txt`](requirements.txt) to easily build the `env`

Following these commands

```
conda create -n env_name python=3.8 # install the base environment
conda activate env_name             # activate the environment
pip install -r requirements.txt     # install the environment
```

## How to setup the program
Download or clone the project into your desired directory (denoted as `/path/to/folder`). Then, add the project directory `cg-dsDNA-openmm` and its tools subdirectory `cg-dsDNA-openmm/tools` to your environment PATH. In addition, ensure that the files `main.py` and `tools/cleanPDB.py` are marked as executable.

```
export PATH=/path/to/folder/cg-dsDNA-openmm/:$PATH
export PATH=/path/to/folder/cg-dsDNA-openmm/tools/:$PATH
```




## Cookbook

### Step 1: Structure (Monomer) cleanup
The structure should be pre-processed before converting into coarse-grained resolution, we suggest to use `pdb4amber`, a command provided by AmberTools.

The main purpose here is to retain the chain ID and make the residue ID consecutive, both two identifier are necessary for the following structural construction. Especially for coarse graining double strand DNA (dsDNA).

If we did not specify `--noter`, the dsDNA (whole species) will be splited into two parts with `TER` added between consecutive chains (that is not what we expected)
```
pdb4amber -i model.pdb -o reorder.pdb --noter
```


### Step 2: Complex construction
At this stage, you may resort to external programs, such as `packmol` or `gromacs`. Here, we detailed the procedure using these two programs.

***Using `packmol`***

```
add_amber_ter # this command is necessary for the following coarse-graining

tolerance 2.0
filetype pdb
output mixture.pdb

structure species_1.pdb  # specify the location of composition 1
  number 100             # the number you want to put in the box
  inside box 0. 0. 0. 1000. 1000. 1000.  # box info
end structure

structure species_2.pdb  # specify the location of composition 2
  number 25
  inside box 0. 0. 0. 1000. 1000. 1000. 
end structure
```

If `packmol` successfully generates the complex, the resulting structure can be directly used in the next stage without further modification. 

(***This approach is efficient for small systems, but not for much larger ones.***)


***Using `gmx insert-molecules`***

``` 
gmx insert-molecules -ci species_1 -o complex1.pdb -nmol 100 -box 100 100 100 # insert species 1 into specified box
gmx insert-molecules -ci species_2 -f complex1.pdb -o complex1_1.pdb -nmol 25 -box 100 100 100 # insert species 2 into specified box

# you can insert other species following the previous command
```

Once done with `gromacs`, we need to insert `TER` at the end of each molecule.

For your convenience, we provide a script to add `TER` automatically.
```
cleanPDB.py -i md.conf # we will introduce how to prepare md.conf for such simulation
```

### Step 3: configure preparation
Currently, a template is provided and stored in the `template/md.conf` subdirectory.

We will explain this configure line by line.

```
# system info
protein = 10
DNA     = 10
fn_protein =  ../test/abeta.pdb ; file for the monomer protein
fn_DNA     =  ../test/JCTC_dsDNA.pdb ; file for the monomer DNA

# simulation parameters
xyz = ../test/Yi.pdb ; file for the assembled systems
interval = 10_000
timestep = 0.01  ; in the unit of ps
temperature = 300 ; in the unit of K
dielectric  = 80  ; useful when controlling 
				  ;electrostatic interactions
cellVector1 = 1000,0,0
cellVector2 = 0,1000,0
cellVector3 = 0,0,1000
cellOrigin  = 0,0,0

# md control
md          = md ; options include `sp` and `rst`
				; sp for single point estimation
				; rst for restarting md simulations (on going)

# machine properties
platform = CPU
DeviceIndex=0,1
Precision=mixed


# output option
trajfile = trajtest.dcd   ; trajectory file for save trajectories
rstfile  = checkpoint.chk ; checkpoint file for restarting the simulation
log = log.openmm
```

### Step 4: Running simulations

```
main.py -i md.conf -v -lammps 
```

- `-h` or `--help`: show help message and exit
- `-i` or `--i_conf`: feed in the `md.conf` prepared in `Step 3`
- `-v` or `--verbose`: print out more information
- `-lammps` or `--dumplammps`: output `in.data` and `lammps.in` for running lammps
