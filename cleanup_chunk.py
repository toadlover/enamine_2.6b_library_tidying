#this script serves to clean a single chunk (0-530)
#it removes all extraneous shapedb data that is currently saved in the folder
#this then goes through and looks at the ligands in the split_new_named.sdf file and determines if ligands contain unwanted functional groups/atoms
#unwanted atoms/functional groups: esters, nitriles, iodine, bromine
#this also writes a blacklist file to note what ligands were removed and why/which functional groups they had
#the blacklist file will exist for the chunk and list all ignored ligands from all subchunks within the chunk. It is needed because we have no effective way to remove ligands from the shapedb data that is existing, and we do not want to recreate the shapedb data

#imports
import os,sys
from rdkit import Chem

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
os.chdir("/pi/summer.thyme-umw/ari_enamine_conformer_library_tidying/conformer_library_space/" + working_chunk)

#save the working chunk location
working_chunk_location = os.getcwd()

#begin writing a chunk blacklist file
blacklist_file = open("blacklist_file.csv", "w")

#write a header
blacklist_file.write("ligand,chunk,subchunk,splitfile,reasons,SMILES\n")

#iterate over every subchunk folder in the directory to perform the operations on
for r,d,f in os.walk(working_chunk_location):
	for dire in d:
		if r == working_chunk_location:
			#work on directory if it is located at the top of the chunk directory
			print("Processing: " + str(dire))