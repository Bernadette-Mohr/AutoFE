#!/bin/bash

# email on start, end, and abortion
#SBATCH --mail-type=ALL
#SBATCH --mail-user= "================> your email <================"
# echo "Job $SLURM_JOBID started at `date`" | mail $USER -s "Job $SLURM_JOBID"

#SBATCH -t 90:00:00

#SBATCH --job-name=mol-NR
#SBATCH --output=mol-NR.out
#SBATCH --error=mol-NR.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH -p normal


module purge
module load "================> appropriate module <================"
module load "================> GROMACS version 2018.6 <================"

# in case anaconda is used to set python version and modules. Location of conda environment may be different.
source $HOME/.bashrc
conda activate py37

NTASKS=$((1 * $SLURM_JOB_NUM_NODES))
echo "ntasks: "$NTASKS

# SET NUMBER OF OPENMP THREADS
OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "omp-num-threads: "$OMP_NUM_THREADS

# GET NODE JOB WAS ALLOCATED TO
echo "node: "$SLURM_NODELIST

# SET NUMBER OF LAMBDA STEPS
STEPS=$((N_STEPS - 1))

mkdir -p "$TMPDIR"/molecule_NR           #create a temporary directory for my job
mkdir -p "$TMPDIR"/molecule_NR/scripts

cp $HOME/AutoFE/analysis/analysis.py $HOME/AutoFE/analysis/position_probabilities.py $HOME/AutoFE/analysis/process_FEs.py "$TMPDIR"/molecule_NR/scripts/

for env in CDL2 POPG OCTANE WATER; do
    mkdir -p "$TMPDIR"/molecule_NR/$env
    cp -r $env/molecule_NR/* "$TMPDIR"/molecule_NR/$env   #copy input files (eg *.tpr, *.pdb) to scratch memory
done

if [ ! -f "$TMPDIR"/molecule_NR/OCTANE/OCTANE-molecule_NR-39/prod-39.gro ]; then

    cd "$TMPDIR"/molecule_NR/POPG/

    for i in $(seq 0 $STEPS); do

        j=$(printf "%02d" $i)
        k=$(printf "%02d" $((i-1)))

        cp POPG-molecule_NR-$k/prod-$k.gro POPG-molecule_NR-$j/init-$j.gro
        cd POPG-molecule_NR-$j/

        # production run
        if [ ! -f prod-$j.gro ]; then
            gmx_mpi grompp -f prod-$j.mdp -c init-$j.gro -p system.top -n index.ndx -maxwarn 3 -o prod-$j -po prod-run.mdp
            mpirun -np $NTASKS gmx_mpi mdrun -ntomp $OMP_NUM_THREADS -deffnm prod-$j
        fi

        cd ../

    done

    python ../scripts/position_probabilities.py -env $PWD
    is_interfacial=$?

    if [ $is_interfacial -eq 1 ]; then
        
        python ../scripts/analysis.py -env $PWD -s 3000
        
        cp -r POPG-molecule_*-* *.png *.pickle *.txt $HOME/WORK_DIR/POPG/molecule_NR/ #copy output files from scratch to your working directory
        cp -r "$TMPDIR"/molecule_NR/FEs-molecule_NR.pickle $HOME/WORK_DIR/
        cd $HOME/WORK_DIR/

    else

        cp -r POPG-molecule_*-* *.png *.pickle *.txt $HOME/WORK_DIR/POPG/molecule_NR/ #copy output files from scratch to your working directory
        cd "$TMPDIR"/molecule_NR/WATER/

        for i in $(seq 0 $STEPS); do

            j=$(printf "%02d" $i)
            k=$(printf "%02d" $((i-1)))

            cp WATER-molecule_NR-$k/prod-$k.gro WATER-molecule_NR-$j/init-$j.gro
            cd WATER-molecule_NR-$j/

            # production run
            if [ ! -f prod-$j.gro ]; then
                gmx_mpi grompp -f prod-$j.mdp -c init-$j.gro -p system.top -n index.ndx -maxwarn 3 -o prod-$j -po prod-run.mdp
                mpirun -np $NTASKS gmx_mpi mdrun -ntomp $OMP_NUM_THREADS -deffnm prod-$j
            fi
        
            cd ../

        done

        cp -r WATER-molecule_*-* $HOME/WORK_DIR/WATER/molecule_NR/ #copy output files from scratch to your working directory
    
        cd "$TMPDIR"/molecule_NR/OCTANE/

        for i in $(seq 0 39); do
    
            j=$(printf "%02d" $i)
            k=$(printf "%02d" $((i-1)))

            cp OCTANE-molecule_NR-$k/prod-$k.gro OCTANE-molecule_NR-$j/init-$j.gro
            cd OCTANE-molecule_NR-$j/

            # production run
            if [ ! -f prod-$j.gro ]; then
                gmx_mpi grompp -f prod-$j.mdp -c init-$j.gro -p system.top -n index.ndx -maxwarn 3 -o prod-$j -po prod-run.mdp
                mpirun -np $NTASKS gmx_mpi mdrun -ntomp $OMP_NUM_THREADS -deffnm prod-$j
            fi

            cd ../

        done
    
        cp -r OCTANE-molecule_*-* $HOME/WORK_DIR/OCTANE/molecule_NR/ #copy output files from scratch to your working directory
    
        cd ../

        python "$TMPDIR"/molecule_NR/scripts/analysis.py -env "$TMPDIR"/molecule_NR/POPG "$TMPDIR"/molecule_NR/OCTANE "$TMPDIR"/molecule_NR/WATER -s 3000
    
        python "$TMPDIR"/molecule_NR/scripts/process_FEs.py -a "$TMPDIR"/molecule_NR/FEs-molecule_NR.pickle
        
        cp "$TMPDIR"/molecule_NR/*.pickle $HOME/WORK_DIR/
        cp "$TMPDIR"/molecule_NR/POPG/*.txt $HOME/WORK_DIR/POPG/molecule_NR/
        cp "$TMPDIR"/molecule_NR/WATER/*.txt $HOME/WORK_DIR/WATER/molecule_NR/
        cp "$TMPDIR"/molecule_NR/OCTANE/*.txt $HOME/WORK_DIR/OCTANE/molecule_NR/
        
    fi

    cd $HOME/WORK_DIR/
    echo "end of first part of screening job"    

else

    cp $HOME/WORK_DIR/FEs-molecule_NR.pickle "$TMPDIR"/molecule_NR/
    cd "$TMPDIR"/molecule_NR/CDL2/
    
    for i in $(seq 0 $STEPS); do
    
        j=$(printf "%02d" $i)
        k=$(printf "%02d" $((i-1)))
        
        cp CDL2-molecule_NR-$k/prod-$k.gro CDL2-molecule_NR-$j/init-$j.gro
        cd CDL2-molecule_NR-$j/

        # production run
        if [ ! -f prod-$j.gro ]; then
            gmx_mpi grompp -f prod-$j.mdp -c init-$j.gro -p system.top -n index.ndx -maxwarn 3 -o prod-$j -po prod-run.mdp
            mpirun -np $NTASKS gmx_mpi mdrun -ntomp $OMP_NUM_THREADS -deffnm prod-$j
        fi

        cd ../
    
    done

    cp -r CDL2-molecule_*-* $HOME/WORK_DIR/CDL2/molecule_NR/ #copy output files from scratch to your working directory
        
    cd ../
        
    python "$TMPDIR"/molecule_NR/scripts/analysis.py -env "$TMPDIR"/molecule_NR/CDL2 -s 3000 -a "$TMPDIR"/molecule_NR/FEs-molecule_NR.pickle
    
    python "$TMPDIR"/molecule_NR/scripts/process_FEs.py -a "$TMPDIR"/molecule_NR/FEs-molecule_NR.pickle
        
    cp "$TMPDIR"/molecule_NR/*.pickle $HOME/WORK_DIR/
    cp "$TMPDIR"/molecule_NR/CDL2/*.txt $HOME/WORK_DIR/CDL2/molecule_NR/
    
    cd $HOME/WORK_DIR/
    echo "end of screening job"

fi

