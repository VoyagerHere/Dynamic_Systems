# srun -p a100 -c 32 -t 4320 ~/julia-1.10.2/bin/julia -t 32  fr_discr_sin_110.005.jl
# srun -p a100 -c 32 -t 4320 ~/julia-1.10.2/bin/julia -t 32  fr_discr_sin_110.01.jl
# srun -p a100 -c 32 -t 4320 ~/julia-1.10.2/bin/julia -t 32  fr_discr_sin_110.05.jl
# srun -p a100 -c 32 -t 4320 ~/julia-1.10.2/bin/julia -t 32  fr_discr_sin_330.005.jl
# srun -p a100 -c 32 -t 4320 ~/julia-1.10.2/bin/julia -t 32  fr_discr_sin_330.01.jl
# srun -p a100 -c 32 -t 4320 ~/julia-1.10.2/bin/julia -t 32  fr_discr_sin_330.05.jl

srun -p a100 -c 32 -t 4320 ~/julia-1.10.2/bin/julia -t 32  fr_discr_cos_110.005.jl
srun -p a100 -c 32 -t 4320 ~/julia-1.10.2/bin/julia -t 32  fr_discr_cos_110.01.jl
srun -p a100 -c 32 -t 4320 ~/julia-1.10.2/bin/julia -t 32  fr_discr_cos_110.05.jl
srun -p a100 -c 32 -t 4320 ~/julia-1.10.2/bin/julia -t 32  fr_discr_cos_330.005.jl
srun -p a100 -c 32 -t 4320 ~/julia-1.10.2/bin/julia -t 32  fr_discr_cos_330.01.jl
srun -p a100 -c 32 -t 4320 ~/julia-1.10.2/bin/julia -t 32  fr_discr_cos_330.05.jl