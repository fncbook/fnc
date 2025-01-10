using Plots

default(label="",linewidth=2,fontfamily="Helvetica")
plt = plot(dpi=400,size=(800,220),label="",layout=(1,2),frame=:none);
s = sqrt(pi)*0.18;
x = range(0-s,sqrt(pi)+s,length=500)
plot!(x,sin.(mod.(x,sqrt(pi)).^2),subplot=1,l=:dash)
x = range(0,sqrt(pi),length=500)
plot!(x,sin.(x.^2),subplot=1)
scatter!([0,sqrt(pi)],[0,0],m=(4,:black),subplot=1)

s = 2pi*0.18;
x = range(0-s,2(pi)+s,length=500)
plot!(x,(@. exp(-sin(x))),subplot=2,l=:dash)
x = range(0,2(pi),length=500)
plot!(x,(@. exp(-sin(x))),subplot=2)
scatter!([0,2(pi)],[1,1],m=(4,:black),subplot=2)

savefig(plt,"periodicfun.pdf")
savefig(plt,"periodicfun.svg")