import os,sys
import fnmatch

def directory_sizes_with_pattern(root, pattern):
    logical_total = 0
    disk_total = 0
    seen_inodes = set()

    for dirpath, dirnames, filenames in os.walk(root):
        for name in filenames:
            if not fnmatch.fnmatch(name, pattern):
                continue

            path = os.path.join(dirpath, name)

            try:
                st = os.stat(path)
            except FileNotFoundError:
                continue

            # ls-like size
            logical_total += st.st_size

            # du-like size (avoid double counting hard links)
            inode = (st.st_dev, st.st_ino)
            if inode not in seen_inodes:
                seen_inodes.add(inode)
                disk_total += st.st_blocks * 512

    return logical_total, disk_total

#create file to write to
size_breakdown_file = open("size_breakdown.csv", "w")

size_breakdown_file.write("superchunk,ll_size,du_size\n")

#iterate over all numerical directories in the working directory
for i in range(0,531):
	#only perform if the directory exists
	if os.path.isdir(str(i)):
		print("On superchunk: " + str(i))
		logical_total, disk_total = directory_sizes_with_pattern(root = str(i), pattern = "condensed_params*.tar.gz")
		size_breakdown_file.write(str(i) + "," + str(logical_total) + "," + str(disk_total) + "\n")
		print(str(i) + "," + str(logical_total) + "," + str(disk_total) + "\n")
