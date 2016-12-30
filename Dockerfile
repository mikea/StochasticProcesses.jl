FROM skippa/julia

COPY . StochasticProcesses.jl/
WORKDIR StochasticProcesses.jl/

# RUN ls -al test/runtests.jl
RUN julia -e 'Pkg.clone(pwd()); Pkg.build("StochasticProcesses"); Pkg.test("StochasticProcesses");'

# RUN julia -e ' Pkg.test("StochasticProcesses"; coverage=true)'
# RUN julia -e 'include("test/runtests.jl")'