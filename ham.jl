ENV["MPLBACKEND"]="tkagg";
using QuantumOptics, PyPlot

rc("font", family="serif")

pygui(true)



function Ham(h::Float64,j::Int)
    N=2*j;
    b=SpinBasis(N);
    #vi=normalize(spinup(b) + spindown(b))
    #vi=normalize(Ket(b,fill(1,2*N+1)) + Ket(b,fill(-1,2*N+1)))
    hamiltonian=-(1/N) * (sigmax(b)^2) + h*(sigmaz(b)+(N/2)*identityoperator(b));
    return hamiltonian
end

function echo(hi,hf,t)
    ei,vi=eigenstates(Ham(hi,j),info=false)
    echo=abs(expect(exp(dense(-1im*t*Ham(hf,j))),vi[1]))^2
    return echo
end


j=10
hi=0.5
hf=0.75
t=range(0,40,100)
plot(t, echo.(hi,hf,t))
PyPlot.show()
