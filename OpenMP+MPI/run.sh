#!/bin/bash

# ========== CONFIGURATION ==========
SOURCE_FILE="main.cpp"
EXECUTABLE="out"
REMOTE_PATH="/home/saeed/Desktop/PDC-Project/OpenMP+MPI"
HOSTFILE="hosts.txt"
OMP_THREADS=2
NUM_PROCESSES=4
SLAVE_NODES=("192.168.10.3")  # Add more slave IPs if needed

# ========== STEP 1: COMPILE ==========
echo "[Compiling MPI + OpenMP Code...]"
mpic++ -fopenmp -o $EXECUTABLE $SOURCE_FILE

if [ $? -ne 0 ]; then
  echo "Compilation failed. Exiting."
  exit 1
fi

# ========== STEP 2: COPY EXECUTABLE TO SLAVES ==========
echo "[Copying executable to slave nodes...]"
for NODE in "${SLAVE_NODES[@]}"; do
  scp $EXECUTABLE saeed@$NODE:$REMOTE_PATH
done

# ========== STEP 3: EXPORT OMP ENV ==========
export OMP_NUM_THREADS=$OMP_THREADS
echo "[Set OMP_NUM_THREADS=$OMP_THREADS]"

# ========== STEP 4: RUN WITH MPIRUN ==========
echo "[Running the program with mpirun...]"
mpirun -np $NUM_PROCESSES --hostfile $HOSTFILE $REMOTE_PATH/$EXECUTABLE

