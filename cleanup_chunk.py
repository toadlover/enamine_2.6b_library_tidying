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

#declare patterns of interest to filter out as smarts
patterns = {
    "ester": Chem.MolFromSmarts("C(=O)O[#6]"),
    "nitrile": Chem.MolFromSmarts("C#N"),
    "heavy_halogen": Chem.MolFromSmarts("[Br,I]"),
}

#create a list to hold the name of ligands that are added to the blacklist
blacklist_ligand_names = []

#iterate over every subchunk folder in the directory to perform the operations on
for r,d,f in os.walk(working_chunk_location):
	for dire in d:
		if r == working_chunk_location:
			#work on directory if it is located at the top of the chunk directory
			print("Processing: " + str(dire))

			#move to the directory
			os.chdir(dire)

			#delete the NN files and data csvs that are not needed for the working library
			os.system("rm *NN*.tar.gz *.csv")

			#iterate over each ligand in each sdf file to determine which ligands to cut and add to the blacklist
			for r2,d2,f2 in os.walk(r + "/" + dire):
				for file in f2:
					if file.endswith(".sdf.tar.gz"):
						#unzip the file
						os.system("tar -xzf " + file)

						#read the file
						supplier = Chem.SDMolSupplier(file.split(".tar.gz")[0], removeHs=False)

						#iterate over the file
						for mol in supplier:
						    if mol is None:
						        continue  # skip invalid molecules
						    
						    mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else "Unknown"

						    #test print of name
						    print("Ligand: " + mol_name + " " + Chem.MolToSmiles(mol))

						    #check if the molecule has any blacklist fragments
						    #if so, add the ligand to the blacklist and why

						    #create a blacklist string that notes all labels that the ligand may have, starting with a blank string
						    #make the string
						    blacklist_string = ""

						    for label, smarts in patterns.items():
						        contains_fragment = mol.HasSubstructMatch(smarts)

						        #if this has the fragment, add the type to the blacklist string
						        if contains_fragment:
						        	#add to the blacklist string based on if the string was previously empty of not to note multiple fragments if possible
						        	if blacklist_string == "":
						        		blacklist_string = label
						        	else:
						        		#delimit multiple entries with a semicolor
						        		blacklist_string = blacklist_string + ";" + label

						    #if we are dealing with a blacklisted ligand, add it to the blacklist and to the file
						    if blacklist_string != "":
						    	blacklist_ligand_names.append(mol_name)
						    	#blacklist_file.write("ligand,chunk,subchunk,splitfile,reasons,SMILES\n")
						    	blacklist_file.write(mol_name + "," + working_chunk + "," + dire + "," + file + "," + blacklist_string + "," + Chem.MolToSmiles(mol) + "\n")

						#delete the decompressed file
						os.system("rm " + file.split(".tar.gz")[0])

			#iterate over all condensed_params_and_db_*.tar.gz files to remove blacklisted ligands from the library
			for r2,d2,f2 in os.walk(r + "/" + dire):
				for file in f2:
					if file.startswith("condensed_params_and_db_") and file.endswith(".tar.gz"):

						#decompress the file
						os.system("tar -xzf " + file)

						#derive the subchunk based on the file
						subchunk = file.split(".tar")[0].split("_")[len(file.split(".tar")[0].split("_")) - 1]

						#open a write stream of the *_lig_name_list.txt file to write all ligand conformers that should be retained (i.e. remove blacklist)
						write_lig_name_list = open(r2 + "/" + file.split(".tar.gz")[0] + "/temp_" + dire + "_" + subchunk + "lig_name_list.txt", "w")
						read_lig_name_list = open(r2 + "/" + file.split(".tar.gz")[0] + "/" + dire + "_" + subchunk + "lig_name_list.txt", "r")

						#read the *_lig_name_list.txt file and determine which ligands need to be removed
						for line in read_lig_name_list.readlines():
							mylig = line.split("_")[0]

							#if the ligand is not in the blacklist ligand names, write it to the ligand name list and do not attempt to delete anything
							if mylig not in blacklist_ligand_names:
								write_lig_name_list.write(line)








			#return to the top at end
			os.chdir(working_chunk_location)