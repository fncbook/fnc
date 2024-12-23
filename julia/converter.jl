dir = "/Users/driscoll/Documents/GitHub/fnc/chapter11"
files = [
    "blackscholes.md",
    "methodlines.md",
    "absstab-diffusion.md",
    "stiffness.md",
    "boundaries.md",
]
grabbed = []

for file in files
    println("starting $file")
    lines = readlines(joinpath(dir, file))
    n = findfirst(startswith("###DEMO"), lines)
    while !isnothing(n)
        println("found at line $n")
        stop = findnext(startswith("###END"), lines, n + 1)
        name = match(r"\(demo-([^\)]+)\)=", lines[n+1])[1]
        append!(
            grabbed,
            [
                "(demo-$name-julia)="
                "``````{dropdown} @demo-$name"
                lines[n+4:stop-1]
                "``````"
                ""
            ],
        )
        s = split(
            """
            (demo-$name)=
            ::::{prf:example}
            `````{tab-set}
            ````{tab-item} Julia
            :sync: julia
            :::{embed} #demo-$name-julia
            :::
            ````

            ````{tab-item} MATLAB
            :sync: matlab
            :::{embed} #demo-$name-matlab
            :::
            ````

            ````{tab-item} Python
            :sync: python
            :::{embed} #demo-$name-python
            :::
            ````
            `````
            ::::
            """,
            '\n',
        )
        lines = [lines[1:n-1]; s; lines[stop+1:end]]
        n = findnext(startswith("###DEMO"), lines, n + 1)
    end

    open(file, "w") do file
        return foreach(s -> write(file, s * "\n"), lines)
    end
    println("done with $file")
end

open(joinpath(dir, "grabbed.md"), "w") do file
    return foreach(s -> write(file, s * "\n"), grabbed)
end

##
