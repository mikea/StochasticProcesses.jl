#!/bin/bash -eu
./test.sh

echo "Performance Suite..."
docker run -ti -v $PWD:/workspace skippa/julia -e 'push!(LOAD_PATH, "$(pwd())/src"); include("perf/runperf.jl");' | tee perf/runperf.log
docker run -ti -v $PWD:/workspace skippa/julia --track-allocation=user -e 'push!(LOAD_PATH, "$(pwd())/src"); include("perf/cumsim_memprof.jl");' 
