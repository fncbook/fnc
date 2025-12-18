using Plots,LaTeXStrings
using LinearAlgebra

pgfplotsx()
default(legend=nothing,aspect_ratio=1,linewidth=2)

push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepackage{fontspec}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\setmainfont{TeX Gyre Schola}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepackage{unicode-math}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\setmathfont{texgyreschola-math.otf}")

plt = plot(size=(450,450),aspect_ratio=1,showaxis=false,grid=false,framestyle=:none,xlims=(-4,4),ylims=(-4,4));

A=[-2;1];
b=[3,2];

plot!([-3*A[1],3*A[1]],[-3,3]*A[2],l=(:dash,:black));
plot!([0,b[1]],[0,b[2]],l=:darkred,arrow=true);
x = A\b;y = A*x;
plot!([0,y[1]],[0,y[2]],l=:darkblue,arrow=true);
plot!([y[1], b[1]],[y[2], b[2]],l=:darkblue,arrow=true);

scatter!([0],[0],m=(5,:black,:black));

sq = [0.4im, -0.4+0.4im, -0.4];
sq = y[1]+1im*y[2] .+ sq*sign(y[1]+1im*y[2]);
plot!(real(sq),imag(sq),l=(1.5,:darkgray));

# annotate!(-2,0.25,text("range",:right));
annotate!(-1.7,0.5,text(L"\text{range}(\mathbf{A})",:right));
annotate!(.7,0.9,text(L"\mathbf{b}"));
annotate!(-.1,-.1,text(L"\mathbf{0}",:right,:top));
annotate!(1.,-0.7,text(L"\mathbf{A} \mathbf{x}",:right));
annotate!(2.4,0.4,text(L"\mathbf{b} - \mathbf{A} \mathbf{x}",:left));

plt
savefig(plt,"normaleqns2d.svg")
savefig(plt,"normaleqns2d.pdf")