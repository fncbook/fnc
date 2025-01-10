
function grid_to_tip(lines)
    while true
        start = findfirst(contains("::::{grid}"), lines)
        isnothing(start) && break
        card  = findnext(contains(":::{card}"), lines, start + 1)
        cardstop = findnext(contains(":::"), lines, card + 1)
        stop  = findnext(startswith("::::"), lines, cardstop + 1)
        tip = filter(x -> !all(isspace, x), lines[card+1:cardstop-1])
        splice!(lines, card:stop, ["```{tip}"; ":class: dropdown"; tip ; "```"])
        deleteat!(lines, start)
    end
    return lines
end

# init_cell = ["```{code-cell}", ":tags: [remove-cell]", "include(\"../../../julia/FNC_init.jl\")", "```"]
init_cell = ["```{code-cell}", ":tags: [remove-cell]", "exec(open(\"../../../python/FNC_init.py\").read())", "```"]

lang = "python"
yaml = open("yaml.txt", "w")
for chap in [12]
    println("chapter $chap")
    local stop, start, name, file, output
    file = "$lang/chapter$chap.md"
    dest = "chapter$chap"

    lines = readlines(file)
    header = lines[1:findnext(startswith("---"), lines, 2)]
    start = length(header)
    for section in 1:12
        start = findfirst(startswith("### $chap.$section"), lines)
        isnothing(start) && continue
        println("    section $section")
        mkpath(joinpath(dest, "section$section", lang))
        stop = something(findnext(startswith("### $chap."), lines, start+1), length(lines) + 1)
        while true
            start = something(findnext(startswith("(demo-"), lines, start+1), length(lines) + 1)
            if start >= stop
                break
            end
            fullname = match(Regex("(demo-.*-[^-]*)-$lang"), lines[start])[1]
            name = match(Regex("demo-.*-([^-]*)-$lang"), lines[start])[1]
            if isnothing(name)
                println("Error: $chap.$section, line $start")
                break
            end
            start = findnext(contains("{dropdown}"), lines, start+1)
            output = joinpath(dest, "section$section", lang, "$name.md")
            open(output, "w") do io
                tofile = s -> println(io, s)
                foreach(tofile, header)
                foreach(tofile, init_cell)
                tofile("[**Demo %s**](#$fullname)")
                tofile("")
                finish = findnext(contains("``````"), lines, start+1)
                if isnothing(finish)
                    println("Error: $chap.$section, line $start, $name")
                else
                    text = grid_to_tip(copy(lines[start+1:finish-1]))
                    foreach(tofile, text)
                end
            end
            try
                # run(`jupytext --to notebook --execute $output`)
                println(yaml, "- file: chapter$chap/section$section/$lang/$name.ipynb")
            catch
                println("Error executing: $chap.$section, $name")
            end
        end
    end
end
close(yaml)

##
# for (root, dirs, files) in walkdir("julia/notebooks/chapter13")
#     for file in filter(endswith(".md"), files)
#         println(joinpath(root, file))
#         lines = readlines(joinpath(root, file))
#         open(joinpath(root, file), "w") do io
#             foreach(s -> println(io, s), lines)
#         end
#     end
# end
