using Plots,LaTeXStrings
using LinearAlgebra,FundamentalsNumericalComputation

gr()
default(legend=nothing,linewidth=2)
plt = plot(aspect_ratio=1,showaxis=false,grid=false,xaxis=:none,yaxis=false)

A=[-2;1];
b=[3,2];

plot!([-3*A[1],3*A[1]],[-3,3]*A[2],l=(:dash,:black))
plot!([0,b[1]],[0,b[2]],l=:darkred)
x = A\b;y = A*x;
plot!([0,y[1]],[0,y[2]],l=:darkblue)
plot!([b[1], y[1]],[b[2], y[2]],l=:darkblue)

scatter!([0],[0],m=(5,:black,:black))

sq = [0.4im, -0.4+0.4im, -0.4];
sq = y[1]+1im*y[2] .+ sq*sign(y[1]+1im*y[2]);
plot!(real(sq),imag(sq),l=:darkgray)
xlims!(-3.5,3.5)
ylims!(-3.5,3.5)

annotate!(-2,0.25,text("range",:right))
annotate!(-1.93,0.25,text(L"(\!\mathbf{A})",:left))
annotate!(.7,0.9,text(L"\mathbf{b}"))
annotate!(1.,-0.8,text(L"\mathbf{A} \mathbf{x}",:right))
annotate!(2.4,0.4,text(L"\mathbf{b} - \mathbf{A} \mathbf{x}",:left))


savefig(plt,"normaleqns2d.svg")