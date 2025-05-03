using CairoMakie, LaTeXStrings, Colors

a, b = -3, 1
c, d = -2.1, 1.9
m, n = 3, 5

x = range(a, b, length=m+1)
y = range(c, d, length=n+1)
grid = [Point2f(x, y) for x in x, y in y]

fig = Figure(size=(600, 400))
ax = Axis(fig[1, 1])

for i in 1:m+1
    text!(x[i], -0.1, text=latexstring("x_$(i-1)"), fontsize=20, align=(:center, :top))
    lines!([x[i], x[i]], [c, d], color=:lightgray, linewidth=0.75)
end

for j in 1:n+1
    text!(0.1, y[j], text=latexstring("y_$(j-1)"), fontsize=20, align=(:left, :center))
    lines!([a, b], [y[j], y[j]], color=:lightgray, linewidth=0.75)
end


scatter!(vec(grid), markersize=12, color=Gray(0.3))

hlines!([0], linewidth=3, color=:black)
vlines!([0], linewidth=3, color=:black)

scatter!(x, 0*x, markersize=15, color=:blue)
scatter!(0*y, y, markersize=15, color=:red)

text!(x[2]+0.06, y[5]+0.05, text=L"(x_1, y_4)", fontsize=20)


hidespines!(ax)
hidedecorations!(ax)
ax.xticks = x
ax.yticks = y
# ax.xgridvisible = ax.ygridvisible = true


fig
