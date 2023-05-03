ENV["MPLBACKEND"]="tkagg";
using QuantumOptics, PyPlot

rc("font", family="serif")

pygui(true)



function Ham(α::Float64,j::Int)
    N=2*j;
    b=SpinBasis(N);
    hamiltonian=-(1/N) * (sigmax(b)^2) + α*(sigmaz(b)+(N/2)*identityoperator(b));
    e=eigenenergies(dense(hamiltonian))/N
end

j=5
α=range(0,2,100)
plot(α, Ham.(α,j))
PyPlot.show()