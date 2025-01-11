using Plots,LaTeXStrings
using LinearAlgebra

pgfplotsx()
default(legend=nothing,aspect_ratio=1,linewidth=2)

push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepackage{fontspec}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\setmainfont{TeX Gyre Schola}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepackage{unicode-math}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\setmathfont{texgyreschola-math.otf}")

plt = plot(size=(600,300),showaxis=false,layout=(1,2),framestyle=:none);

u = .3 - .2im;
v = exp(0.1im)*u;

plot!(xlims=(-0.1,.45),ylims=(-.25,.3),aspect_ratio=1,subplot=1);
plot!([0,real(u)],0.1.+[0,imag(u)],l=:darkred,arrow=true,subplot=1);
plot!([0, real(v)],0.1.+[0, imag(v)],l=:darkblue,arrow=true,subplot=1);
plot!(real([v, u]),0.1.+imag([v, u]),l=:black,arrow=true,subplot=1);

u = .3 - .2im;
v = 1im*u;
plot!(xlims=(-0.1,.45),ylims=(-.25,.3),aspect_ratio=1,subplot=2);
plot!([0, real(u)],[0, imag(u)],l=:darkred,arrow=true,subplot=2);
plot!([0 ,real(v)],[0, imag(v)],l=:darkblue,arrow=true,subplot=2);
plot!(real([v, u]),imag([v, u]),l=:black,arrow=true,subplot=2);

annotate!(.15,-0.02,text(L"\mathbf{u}",:right),subplot=1);
annotate!(.18,0.02,text(L"\mathbf{v}",:left),subplot=1);
annotate!(.325,-0.09,text(L"\mathbf{u}-\mathbf{v}",:left),subplot=1);
annotate!(.15,-0.12,text(L"\mathbf{u}",:right),subplot=2);
annotate!(.02,0.15,text(L"\mathbf{v}",:left),subplot=2);
annotate!(.25,0.12,text(L"\mathbf{u}-\mathbf{v}",:left),subplot=2);

plt
savefig("nonorthogonal.svg")
savefig("nonorthogonal.pdf")
