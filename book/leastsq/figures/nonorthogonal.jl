using Plots,LaTeXStrings
using LinearAlgebra,FundamentalsNumericalComputation

gr()
default(legend=nothing,linewidth=2)
plt = plot(aspect_ratio=1,showaxis=false,grid=false,xaxis=:none,yaxis=false,layout=(1,2))

u = .3 - .2im;
v = exp(0.1im)*u;

plot!([0,real(u)],0.1.+[0,imag(u)],l=:darkred,subplot=1)
plot!([0, real(v)],0.1.+[0, imag(v)],l=:darkblue,subplot=1)
plot!(real([u, v]),0.1.+imag([u, v]),l=:black,subplot=1)
xlims!(0 ,.34) 
ylims!(-.25, .3)

u = .3 - .2im;
v = 1im*u;
plot!([0, real(u)],[0, imag(u)],l=:darkred,subplot=2)
plot!([0 ,real(v)],[0, imag(v)],l=:darkblue,subplot=2)
plot!(real([u, v]),imag([u, v]),l=:black,subplot=2)
xlims!(0 ,.34) 
ylims!(-.25, .3)

annotate!(.15,-0.02,text(L"\mathbf{u}",:right),subplot=1)
annotate!(.18,0.02,text(L"\mathbf{v}",:left),subplot=1)
annotate!(.325,-0.09,text(L"\mathbf{u}-\mathbf{v}",:left),subplot=1)
annotate!(.15,-0.12,text(L"\mathbf{u}",:right),subplot=2)
annotate!(.06,0.15,text(L"\mathbf{v}",:left),subplot=2)
annotate!(.27,0.12,text(L"\mathbf{u}-\mathbf{v}",:left),subplot=2)

savefig("nonorthogonal.svg")