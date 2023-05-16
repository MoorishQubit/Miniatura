ENV["MPLBACKEND"]="tkagg"
using QuantumOptics, PyPlot, Arpack, ProgressMeter, LinearAlgebra, NPZ

pygui(true)
rc("font", family="serif")



BLAS.set_num_threads(20)

function Ham(h::Float64,j::Int)
    N=2*j;
    b=SpinBasis(N);
    hamiltonian=-(1/N) * (sigmax(b)^2) + h*(sigmaz(b)+(N/2)*identityoperator(b));
    return hamiltonian
end

function echo(hi,hf,j,t)
    data=[]
    p = Progress(length(t), "Calculating LE for N=$(2*j)")
    ei,vi=eigenstates(Ham(hi,j),nev=1,which=:SM, info=false)
    Threads.@threads for i ∈ t
        l=abs(expect(exp(-1im*dense(Ham(hf,j))*i),vi)[1])^2
        push!(data,l)
        next!(p)
    end
    finish!(p)
    return data
end

j=100; hi=0.5; hf=0.75;
t=range(0,5,50)
# echo(hi,hf,j,t)
PyPlot.plot(t,echo(hi,hf,j,t))
PyPlot.show()