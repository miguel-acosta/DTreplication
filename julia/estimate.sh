#!/bin/bash
#
#SBATCH --account=sscc          # The account name for the job.
#SBATCH --job-name=lDTMHfull     # The job name.
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH -c  4                  # The number of cpu cores to use.
#SBATCH --time=47:59:59        # The time the job will take to run.
#SBATCH --mem-per-cpu=4gb     # The memory the job will use per cpu core.
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jma2241@columbia.edu
#SBATCH -e /rigel/home/jma2241/coursework/uribe/error_HPC.txt
#SBATCH -o /rigel/home/jma2241/coursework/uribe/output_HPC.txt

module load julia/0.6.1
module load anaconda/3-4.2.0
export JULIA_PKGDIR=/rigel/sscc/users/jma2241/.julia/



cd /rigel/home/jma2241/coursework/uribe/DTreplication/julia/
# Arguments: See do file

julia estimation.jl 
