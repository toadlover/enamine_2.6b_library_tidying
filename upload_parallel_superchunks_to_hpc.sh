#!/bin/bash

# Usage:
#   ./parallel_scp.sh START_INDEX END_INDEX PARALLEL_LIMIT
#
# Example:
#   ./parallel_scp.sh 0 530 4

START=$1
END=$2
PARALLEL=$3

# ---- UPDATE THESE ----
KEY="/store/projects/ari/id_rsa"
DEST="ari.ginsparg-umw@hpc.umassmed.edu:/pi/summer.thyme-umw/enamine-REAL-2.6billion"
# -----------------------

if [[ -z "$START" || -z "$END" || -z "$PARALLEL" ]]; then
    echo "Usage: $0 START_INDEX END_INDEX PARALLEL_LIMIT"
    exit 1
fi

# Build list of directory names
DIRS=()
for ((i=START; i<=END; i++)); do
    if [[ -d "$i" ]]; then
        DIRS+=("$i")
    else
        echo "Skipping $i (not a directory)"
    fi
done

export KEY
export DEST

printf "%s\n" "${DIRS[@]}" | \
    xargs -I{} -P "$PARALLEL" bash -c '
        FOLDER="{}"
        echo "Uploading directory: $FOLDER"
        scp -i "$KEY" -dr "$FOLDER" "$DEST"
        STATUS=$?
        if [[ $STATUS -ne 0 ]]; then
            echo "ERROR: SCP failed for $FOLDER"
        else
            echo "Completed: $FOLDER"
        fi
    '
