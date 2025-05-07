#!/bin/bash

# ========== CONFIGURATION ==========
EXECUTABLE="out"
REMOTE_PATH="/home/saeed/Desktop/PDC-Project/MPI"
HOSTFILE="hosts.txt"
NUM_PROCESSES=4
SLAVE_NODES=("192.168.10.3")  # Add more if needed

# ========== STEP 1: COPY EXECUTABLE TO SLAVES ==========
echo "[Copying executable to slave nodes...]"
for NODE in "${SLAVE_NODES[@]}"; do
    scp $EXECUTABLE saeed@$NODE:$REMOTE_PATH/
done

# ========== STEP 2: RUN WITH MPIRUN ==========
echo "[Running MPI program with $NUM_PROCESSES processes...]"
mpirun -np $NUM_PROCESSES --hostfile $HOSTFILE $REMOTE_PATH/$EXECUTABLE

