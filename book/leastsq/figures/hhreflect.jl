using Plots,LaTeXStrings
using LinearAlgebra,FundamentalsNumericalComputation

z = [-1,3];
v = z - norm(z)*[1,0];
P = I - 2*(v*v')/dot(v,v)
Pz = P*z;

gr()
default(legend=nothing,linewidth=1)
plt = plot(aspect_ratio=1,showaxis=false,grid=false,xaxis=:none,yaxis=false)

plot!([0,z[1]],[0,z[2]],aspect_ratio=1)
plot!([Pz[1],z[1]],[Pz[2],z[2]])
plot!([0,Pz[1]],[0,Pz[2]])

plot!([-2,2]*v[2],[2 ,-2]*v[1],l=(:dash,:black))

scatter!([0],[0],m=(5,:black,:black))
xlims!(-1,4)
ylims!(-1, 3)

annotate!(-0.68,1.5,text(L"\mathbf{z}"))
annotate!(2,1,text(L"\mathbf{v}"))
annotate!(1.5,-0.3,text(L"\|\mathbf{z}\|\,\mathbf{e}_1"))
annotate!(-0.2,0,text(L"0"))
annotate!(1.6,1.9,text(L"S"))

savefig(plt,"hhreflect.svg")