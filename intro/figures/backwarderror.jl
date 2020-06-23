using Plots,LaTeXStrings
pgfplotsx()
default(legend=nothing,aspect_ratio=1)

plt = plot(showaxis=false,grid=false,xaxis=:none,yaxis=false)

plot!([0,.4, .4, 0],[.9, .9, .1, .1],l=nothing,seriestype=:shape,fillcolor=RGB(1,1,.7))
plot!([1, .6, .6, 1],[.9, .9, .1, .1],l=nothing,seriestype=:shape,fillcolor=RGB(1,1,.7))
 
plot!([.86, .7],[.75, .2],l=(:dash,2,:black))
plot!([.2, .3],[.75, .5],l=(:dash,2,:black))
 
plot!([.2, .86],[.75,.75],l=(2,:darkblue),m=5,markercolor=:black,arrow=true)
plot!([.2, .7],[.75,.2],l=(2,:darkred),arrow=true)
plot!([.3, .7],[.5,.2],l=(2,:darkblue),m=5,markercolor=:black,arrow=true)

annotate!(.875,.75,text(L"y",:left))
annotate!(.72,.18,text(L"\tilde{y}",:left))
annotate!(.18,.77,text(L"x",:right))
annotate!(.29,.47,text(L"\tilde{x}",:right))
annotate!(.5,.79,text(L"f(x)",:top))
annotate!(.51,.52,text(L"\tilde{f}(x)",:bottom))
annotate!(.48,.25,text(L"f(\tilde{x})",:top))
annotate!(.79,.475,text("error",:left))
annotate!(.26,.6,text("backward error",:right))
annotate!(.2,.95,text("data",:center))
annotate!(.8,.95,text("results",:center))

xlims!(0,1)
ylims!(0,1)

savefig(plt,"backwarderror.svg")