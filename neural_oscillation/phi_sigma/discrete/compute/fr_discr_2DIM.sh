srun -p a100 -c 64 -t 4320 ~/julia-1.10.2/bin/julia -t 64 fr_discr_sin_2DIM_110.005.jl
srun -p a100 -c 64 -t 4320 ~/julia-1.10.2/bin/julia -t 64 fr_discr_sin_2DIM_110.01.jl
srun -p a100 -c 64 -t 4320 ~/julia-1.10.2/bin/julia -t 64 fr_discr_sin_2DIM_110.05.jl
srun -p a100 -c 64 -t 4320 ~/julia-1.10.2/bin/julia -t 64 fr_discr_sin_2DIM_330.005.jl
srun -p a100 -c 64 -t 4320 ~/julia-1.10.2/bin/julia -t 64 fr_discr_sin_2DIM_330.01.jl
srun -p a100 -c 64 -t 4320 ~/julia-1.10.2/bin/julia -t 64 fr_discr_sin_2DIM_330.05.jl

srun -p a100 -c 64 -t 4320 ~/julia-1.10.2/bin/julia -t 64 fr_discr_cos_2DIM_110.005.jl
srun -p a100 -c 64 -t 4320 ~/julia-1.10.2/bin/julia -t 64 fr_discr_cos_2DIM_110.01.jl
srun -p a100 -c 64 -t 4320 ~/julia-1.10.2/bin/julia -t 64 fr_discr_cos_2DIM_110.05.jl
srun -p a100 -c 64 -t 4320 ~/julia-1.10.2/bin/julia -t 64 fr_discr_cos_2DIM_330.005.jl
srun -p a100 -c 64 -t 4320 ~/julia-1.10.2/bin/julia -t 64 fr_discr_cos_2DIM_330.01.jl
srun -p a100 -c 64 -t 4320 ~/julia-1.10.2/bin/julia -t 64 fr_discr_cos_2DIM_330.05.jl