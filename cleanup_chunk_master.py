#this script serves to clean a single chunk (0-530)
#it removes all extraneous shapedb data that is currently saved in the folder
#this then goes through and looks at the ligands in the split_new_named.sdf file and determines if ligands contain unwanted functional groups/atoms
#unwanted atoms/functional groups: esters, nitriles, iodine, bromine
#this also writes a blacklist file to note what ligands were removed and why/which functional groups they had
#the blacklist file will exist for the chunk and list all ignored ligands from all subchunks within the chunk. It is needed because we have no effective way to remove ligands from the shapedb data that is existing, and we do not want to recreate the shapedb data

#imports
import os,sys
from rdkit import Chem

#starting location
starting_location = os.getcwd()

#take in a single value as a chunk folder
working_chunk = sys.argv[1]

#bool to identify if the chunk is valid or not
chunk_id_valid = False

#confirm this points to a chunk, otherwise exist
for i in range(0,531):
	if int(working_chunk) == i:
		chunk_id_valid = True
		break

if chunk_id_valid == False:
	print("Inputted chunk " + str(working_chunk) + " is not valid. Quitting.")
	exit()


#print confirmation
print("Working on chunk: " + str(working_chunk))

#move to the chunk
os.chdir("/pi/summer.thyme-umw/enamine-REAL-2.6billion/" + working_chunk)

#save the working chunk location
working_chunk_location = os.getcwd()



#declare patterns of interest to filter out as smarts
patterns = {
    #not doign blanket ester filter
    #"ester": Chem.MolFromSmarts("C(=O)O[#6]"),
    
    "nitrile": Chem.MolFromSmarts("C#N"),
    "heavy_halogen": Chem.MolFromSmarts("[Br,I]"),


    #ester types that we want to remove as they are generally unstable and liable to fracture and form unwanted byproducts
    "simple_alkyl_ester": Chem.MolFromSmarts("[CX3](=O)[OX2][CH3,CH2]"),
    "benzylic_ester": Chem.MolFromSmarts("C(=O)O[CH2][c]"),
    "hetero_activated_ester": Chem.MolFromSmarts("C(=O)O[CH2][N,O,S]"),
    "nitrile_adjacent_ester": Chem.MolFromSmarts("C(=O)O[CH2]C#N"),
    "nitro_adjacent_ester": Chem.MolFromSmarts("C(=O)O[CH2][N+](=O)[O-]"),
    "heterocycle_ester": Chem.MolFromSmarts("[n,o,s;H0]C(=O)O"),
    "acylal": Chem.MolFromSmarts("C(=O)(O)O"),
}

#ester types that we effectively want to keep
#generally listing the types, but not planning to implement a check for these
#"lactone": Chem.MolFromSmarts("O=C1OCC1")
#"lactone_general": Chem.MolFromSmarts("[$([CX3](=O)O),$([CX3](=O)[O-])]=O;R")
#"aryl_ester": Chem.MolFromSmarts("c-O-C(=O)")
#"tert_ester": Chem.MolFromSmarts("C(=O)O[C;X4]([C])[C]")



#iterate over every subchunk folder in the directory to perform the operations on
for r,d,f in os.walk(working_chunk_location):
	for dire in d:
		if r == working_chunk_location:
			#work on directory if it is located at the top of the chunk directory
			print("Processing: " + str(dire))

			#move to the directory
			os.chdir(dire)

			print("working in: " + str(os.getcwd()))

			#bsub job throttle to make sure we do not exceed our local limit
			#write the length of the bjobs queue to this current location
			os.system("bjobs | wc -l > bjobs_length.txt")
			job_count = 0
			with open("bjobs_length.txt") as f:
				job_count = int(f.read().strip())
			while job_count > 50:
				#sleep for 1 second to not overburden the system
				os.system("sleep 1")
				os.system("bjobs | wc -l > bjobs_length.txt")
				with open("bjobs_length.txt") as f:
					job_count = int(f.read().strip())
			#remove the length file to avoid clutter
			os.system("rm bjobs_length.txt")

			#cut things off here and prepare to run this in a bsub job in parallel for cleanup_chunk_sub.py
			os.system("bsub -q long -W 8:00 \"python /pi/summer.thyme-umw/ari_enamine_conformer_library_tidying/enamine_2.6b_library_tidying/cleanup_chunk_sub.py " + working_chunk + " " +  dire + " \"")

			#return to the top at end
			os.chdir(working_chunk_location)