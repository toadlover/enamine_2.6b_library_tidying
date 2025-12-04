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
dire = sys.argv[2]

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




#work on directory if it is located at the top of the chunk directory
print("Processing: " + str(dire))


print("working in: " + str(os.getcwd()))

#delete the NN files and data csvs that are not needed for the working library
os.system("rm *NN*.tar.gz *.csv")

#begin writing a chunk blacklist file
blacklist_file = open("blacklist_file.csv", "w")

#write a header
blacklist_file.write("ligand,superchunk,chunk,subchunk_splitfile,reasons,SMILES\n")

#create a list to hold the name of ligands that are added to the blacklist
blacklist_ligand_names = {}



#iterate over each ligand in each sdf file to determine which ligands to cut and add to the blacklist
for r2,d2,f2 in os.walk(os.getcwd()):
	for file in f2:
		if file.endswith(".sdf.tar.gz"):
			#unzip the file
			os.system("tar -xzf " + file)

			#derive the subchunk based on the file
			subchunk = file.split(".tar")[0].split("_")[len(file.split(".tar")[0].split("_")) - 1]

			blacklist_ligand_names[subchunk] = []

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
			    	blacklist_ligand_names[subchunk].append(mol_name)
			    	#blacklist_file.write("ligand,chunk,subchunk,splitfile,reasons,SMILES\n")
			    	blacklist_file.write(mol_name + "," + working_chunk + "," + dire + "," + file + "," + blacklist_string + "," + Chem.MolToSmiles(mol) + "\n")

			#delete the decompressed file
			os.system("rm " + file.split(".tar.gz")[0])

print(blacklist_ligand_names)

#iterate over all condensed_params_and_db_*.tar.gz files to remove blacklisted ligands from the library
for r2,d2,f2 in os.walk(os.getcwd()):
	for file in f2:
		if file.startswith("condensed_params_and_db_") and file.endswith(".tar.gz"):

			print(file)

			#decompress the file
			os.system("tar -xzf " + file)

			#delete the split_confs file, we do not need an extra
			os.system("rm " + r2 + "/" + file.split(".tar.gz")[0] + "/split_new_named*sdf")


			#derive the subchunk based on the file
			subchunk = file.split(".tar")[0].split("_")[len(file.split(".tar")[0].split("_")) - 1]

			#open a write stream of the *_lig_name_list.txt file to write all ligand conformers that should be retained (i.e. remove blacklist)
			write_lig_name_list = open(r2 + "/" + file.split(".tar.gz")[0] + "/temp_" + dire + "_" + subchunk + "_lig_name_list.txt", "w")
			read_lig_name_list = open(r2 + "/" + file.split(".tar.gz")[0] + "/" + dire + "_" + subchunk + "_lig_name_list.txt", "r")

			#read the *_lig_name_list.txt file and determine which ligands need to be removed
			for line in read_lig_name_list.readlines():
				mylig = line.split("_")[0]

				#if the ligand is not in the blacklist ligand names, write it to the ligand name list and do not attempt to delete anything
				if mylig not in blacklist_ligand_names[subchunk]:
					write_lig_name_list.write(line)
				else:
					#delete the shorthand params file
					os.system("rm -f " + r2 + "/" + file.split(".tar.gz")[0] + "/single_conf_params/" + mylig + "_*txt")

			#overwrite the new ligand name list with the old
			os.system("mv " + r2 + "/" + file.split(".tar.gz")[0] + "/temp_" + dire + "_" + subchunk + "_lig_name_list.txt " + r2 + "/" + file.split(".tar.gz")[0] + "/" + dire + "_" + subchunk + "_lig_name_list.txt")

			#condense the new trimmed directory
			#enter the location for ease in condensing
			os.system("tar -czf " + file + " " + file.split(".tar.gz")[0])

			#add 5 seconds to sleep
			#os.system("sleep 5")

			#delete the made folder
			#os.system("(sleep 60 && rm -drf " + r2 + "/" + file.split(".tar.gz")[0] + ") &")
			os.system("rm -drf " + r2 + "/" + file.split(".tar.gz")[0])
			#os.system("(sleep 60 && rm -drf " + r2 + "/" + file.split(".tar.gz")[0] + ") ")
			#try to kill in the background 5 minutes later, and if that fails, we will kill it in the end step
			#os.system("(sleep 300 && rm -drf " + r2 + "/" + file.split(".tar.gz")[0] + ") & ")
			#os.system("(sleep 600 && rm -drf " + r2 + "/" + file.split(".tar.gz")[0] + ") & ")
			#os.system("(sleep 60 && rm -drf " + r2 + "/" + file.split(".tar.gz")[0] + ") ")


#end cleanup to remove any potential dangling empty uncompressed condensed_params_and_db_(0-9) directories
for i in range(0,10):

	cur_folder = "condensed_params_and_db_" + str(i)

	print("condensed_params_and_db_" + str(i))

	while os.path.isdir(cur_folder):

		print("attempting to delete")

		os.system("rm -drf " + cur_folder)

		#sleep 2 seconds before trying again or moving on
		os.system("sleep 2")




