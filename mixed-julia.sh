#!/bin/bash

#SBATCH --job-name=opt-syn
#SBATCH --output=out/VRP_%j.out    ## ชื่อไฟล์ Output (%j = Job-ID)
#SBATCH --error=out/VRP_%j.out     ## ชื่อไฟล์ error (%j = Job-ID)
#SBATCH --time=24:00:00          ## เวลาที่ใช้รันงาน
#SBATCH --partition=mixed         ## ระบุ partition ที่ต้องการใช้งาน
#SBATCH --nodes=1               # node count
#SBATCH --ntasks=1              ## จำนวน tasks ที่ต้องการใช้ในการรัน
#SBATCH --cpus-per-task=32      ## จำนวน code ที่ต้องการใช้ในการรัน
#SBATCH --mail-type=END
#SBATCH --mail-user=payakorn.s@cmu.ac.th

module purge                    ## unload module ทั้งหมด เพราะว่าอาจจะมีการ load module ไว้ก่อนหน้านั้น

# source $HOME/.julia/dev/VRPTW/
module load julia
module load gurobi
# grbgetkey ef322eb1-aad9-4953-84e8-97e0490a2ee1
# grbgetkey 38f90f84-fa06-4d6c-9cd0-0cc50f65ddfa
# grbgetkey e74a4d63-7186-4785-85f2-744b9d7447bd /home/payakorn.s/gurobi.lic
# grbgetkey 33a2f816-9e2c-4023-9eff-db03959eb467 /home/payakorn.s/gurobi.lic
# mv /home/payakorn.s/gurobi.lic /home/payakorn.s/.julia/dev/VehicleRoutingSynPre/gurobi.lic

# export GRB_LICENSE_FILE=/home/payakorn.s/.julia/dev/VRPTW.jl/gurobi1.lic
export GRB_LICENSE_FILE=/home/payakorn.s/.julia/dev/VehicleRoutingSynPre/gurobi.lic

# srun python copy_of_atom_10_payakorn.py           ## สั่งรัน code
# srun julia --threads 32 src/script.jl           ## สั่งรัน code
# srun julia --threads 32 src/script2.jl           ## สั่งรัน code
# srun julia --threads 32 --project=/home/payakorn.s/.julia/dev/VRPTW.jl src/script3.jl           ## สั่งรัน code
# srun julia --threads 32 src/script3.jl          ## สั่งรัน code
srun julia --threads 32 src/script2.jl          ## สั่งรัน code
