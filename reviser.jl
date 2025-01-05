chapter = 13
println("chapter $chapter")
for section = 1:12
    dir = "chapter$chapter/section$section"
    file = try
        filter(endswith(".md"), readdir(dir))
    catch
        continue
    end
    isnothing(file) && continue
    println("    section $section")
    target = joinpath(dir, file[1])
    lines = readlines(target)
    start = 0
    while true
        start = findnext(contains("(demo-"), lines, start + 1)
        isnothing(start) && break
        prf = findnext(contains("{prf:example}"), lines, start + 1)
        title = match(r"\{prf:example\}\s(.*)", lines[prf])[1]
        embed = prf
        embed = findnext(contains("{embed}"), lines, embed + 1)
        name = match(r"\{embed\} #demo-(.*)", lines[embed])[1]
        fhyphen = findfirst(isequal('-'), name)
        lhyphen = findlast(isequal('-'), name)
        name = name[fhyphen+1:lhyphen-1]
        newtext = """
    ```````{prf:example}
    ``````{tab-set}
    `````{tab-item} Julia
    :sync: julia
    ````{dropdown} $title
    :open:
    ```{include} julia/$name.ipynb
    ```
    ````
    `````

    `````{tab-item} MATLAB
    :sync: matlab
    ````{dropdown} $title
    :open:
    ```{include} matlab/$name.ipynb
    ```
    ````
    `````

    `````{tab-item} Python
    :sync: python
    ````{dropdown} $title
    :open:
    ```{include} python/$name.ipynb
    ```
    ````
    `````
    ``````
    ```````
        """
        brace = findfirst(isequal('{'), lines[prf])
        stop = findnext(startswith(lines[prf][1:brace-1]), lines, prf + 1)
        splice!(lines, prf:stop, split(newtext, '\n'))
    end
    open(target, "w") do io
        foreach(s -> println(io, s), lines)
    end
    println("    done $target")
end
