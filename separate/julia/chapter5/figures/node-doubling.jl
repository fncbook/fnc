using CairoMakie

fig = Figure(resolution=(480,300))
ax = Axis(fig[1,1])
z = 1;
for n = 1:4
	x = range(0,1,length=2^n+1)
	y = fill(z,2^n+1)
	lines!([0,1],[z,z],color=:black,linewidth=2)
	scatter!(x[1:2:end],y[1:2:end],color=:black,markersize=7)
	scatter!(x[2:2:end],y[2:2:end],color=:red,markersize=7)
	z -= 0.25
end
limits!(-0.05,1.05,0.05,1.15)
ax.xticks = range(0,1,length=5)
ax.xtickformat = ""
hidespines!(ax)
hidedecorations!(ax)
fig
save("node-doubling.pdf",fig)
# save("node-doubling.svg",fig)