#!/bin/bash -eu

echo "Testing..."
docker run -ti -v $PWD:/workspace skippa/julia -e 'push!(LOAD_PATH, "$(pwd())/src"); include("test/runtests.jl");'