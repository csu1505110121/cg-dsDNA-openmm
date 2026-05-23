# system info
protein = 800
DNA     = 20
fn_protein = ../rigid.pdb  ; file for the monomer protein
fn_DNA     =  ../32bp.pdb ; file for the monomer DNA

# simulation parameters
xyz = complex_addTER.pdb ; file for the assembled systems
interval = 100_000
timestep = 0.01  ; in the unit of ps
temperature = 300 ; in the unit of K
dielectric  = 80  ; useful when controlling 
				  ;electrostatic interactions
C_salt      = 1   ; in the unit of nm


cellVector1 = 150,0,0
cellVector2 = 0,150,0
cellVector3 = 0,0,150
cellOrigin  = 0,0,0
step     = 50_000_000

md = restart
# vecScale = 20

setrigid = True



# ensemble selection and control
ensemble = NVT
ref_bar  = 1.0
ref_bar_freq = 25


# machine properties
platform = CUDA
# Threads = 8
DeviceIndex=0,1
# Precision=double


# output option
rstfile_p = checkpoint.02.chk

monitorfile = monitor.03.txt
trajfile = trajtest.03.dcd   ; trajectory file for save trajectories
rstfile  = checkpoint.03.chk ; checkpoint file for restarting the simulation
log = log.openmm.03

