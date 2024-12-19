chap = "8"

dir = "/Users/driscoll/Documents/GitHub/fnc/chapter$chap"
code_files = [
    "/Users/driscoll/Documents/GitHub/fnc/julia/chapter$chap.md",
    "/Users/driscoll/Documents/GitHub/fnc/matlab/chapter$chap.md",
    "/Users/driscoll/Documents/GitHub/fnc/python/chapter$chap.md"
]
codes = readlines.(code_files);

for file in filter(endswith(".md"), readdir(dir, join=true))
    println(file)
    lines = readlines(file);
    n = findfirst(startswith("(demo-"), lines)
    while !isnothing(n)
        name = match(r"\(demo-([^\)]+)\)=", lines[n])[1]
        # get title from the julia code file
        locs = findfirst.(startswith("(demo-$name"), codes)
        title = match(r".*\{dropdown}\s+(.*)", codes[1][locs[1] + 1])[1]
        lines[n+1] *= " $title"
        for (code, m) in zip(codes, locs)
            code[m+1] = replace(code[m+1], title => "@demo-$name")
        end
        n = findnext(startswith("(demo-"), lines, n+1)
    end

    open(file,"w") do file
        foreach(s->write(file, s*"\n"), lines)
    end
end

for (file, code) in zip(code_files, codes)
    open(file,"w") do file
        foreach(s->write(file, s*"\n"), code)
    end
end
