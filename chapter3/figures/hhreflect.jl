using Plots,LaTeXStrings
using LinearAlgebra

pgfplotsx()
default(legend=nothing,aspect_ratio=1,linewidth=2)

push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepackage{fontspec}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\setmainfont{TeX Gyre Schola}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepackage{unicode-math}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\setmathfont{texgyreschola-math.otf}")

default(legend=nothing,linewidth=2)

##

plt = plot(size=(500,400),aspect_ratio=1,showaxis=false,grid=false,xaxis=((-4,6)),yaxis=((-5,3)),frame=:none);

P = [2,-1]
plot!(P[1]*[-5,5],P[2]*[-5,5],l=(1,:dash,:black),aspect_ratio=1);

v = 0.9*[-1,-2]
x1 = -0.8P
x2 = x1 + v
plot!([x1[1],x2[1]],[x1[2],x2[2]],arrow=true,color=:darkred);
t = x1 + 0.6v + [-0.15,0.05];
annotate!(t[1],t[2],text(L"\mathbf{v}",:right,:bottom));

x = 1.25*[2,-3];
O = 0P
x1 = O
x2 = x1 + x
plot!([x1[1],x2[1]],[x1[2],x2[2]],arrow=true,color=:darkblue);
t = x1 + 0.66(x2-x1) - [0.25,0];
annotate!(t[1],t[2],text(L"\mathbf{x}",:right,:top));

c = dot(x,v) / dot(v,v)
x1 = x2
x2 = x1 - c*v
plot!([x1[1],x2[1]],[x1[2],x2[2]],arrow=true,color=:black);
t = x1 + 0.4*(x2-x1) + [0.15,0];
annotate!(t[1],t[2],text(L"-(\mathbf{v}^T\mathbf{x})\mathbf{v}",:left));

x1 = x2
x2 = x1 - c*v
plot!([x1[1],x2[1]],[x1[2],x2[2]],arrow=true,color=:black);

x1 = O;
plot!([x1[1],x2[1]],[x1[2],x2[2]],arrow=true,color=:darkblue);
t = x1 + 0.66(x2-x1) + [0,0.32];
annotate!(t[1],t[2],text(L"\mathbf{P}\mathbf{x}",:right,:bottom));

scatter!([O[1]],[O[2]],m=(5,:black,:black));
annotate!(O[1]-0.2,O[2]-0.2,text(L"\mathbf{0}",:top,:right));
# annotate!(1.5,-0.3,text(L"\|\mathbf{z}\|\,\mathbf{e}_1"));
# annotate!(-0.2,0,text(L"0"));
# annotate!(1.6,1.9,text(L"S"));

xlims!((-4,6))
ylims!((-5,3))
plt
savefig(plt,"hhreflect.svg")
savefig(plt,"hhreflect.pdf")

