# system info
protein = 800
DNA     = 20
fn_protein = ../rigid.pdb  ; file for the monomer protein
fn_DNA     =  ../32bp.pdb ; file for the monomer DNA

# simulation parameters
xyz = complex.pdb ; file for the assembled systems
interval = 100_000
timestep = 0.01  ; in the unit of ps
temperature = 300 ; in the unit of K
dielectric  = 80  ; useful when controlling 
				  ;electrostatic interactions
cellVector1 = 150,0,0
cellVector2 = 0,150,0
cellVector3 = 0,0,150
cellOrigin  = 0,0,0
step     = 100_000_000

md = md
#vecScale = 20

setrigid = True


# ensemble selection and control
ensemble = NPT
ref_bar  = 1.0
ref_bar_freq = 25



# machine properties
platform = CPU
Threads = 8
#DeviceIndex=0,1
#Precision=mixed


# output option
#rstfile_p = checkpoint.chk

monitorfile = monitor.01.txt
trajfile = trajtest.01.dcd   ; trajectory file for save trajectories
rstfile  = checkpoint.01.chk ; checkpoint file for restarting the simulation
log = log.openmm.01

