#!/bin/bash -eu

echo "Testing..."
docker run -ti -v $PWD:/workspace skippa/julia -e 'push!(LOAD_PATH, "$(pwd())/src"); include("test/runtests.jl");'

echo "Performance Suite..."
docker run -ti -v $PWD:/workspace skippa/julia -e 'push!(LOAD_PATH, "$(pwd())/src"); include("perf/runperf.jl");' | tee perf/runperf.log