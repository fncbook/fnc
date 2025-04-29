chap = 7
dir = "/Users/driscoll/Documents/GitHub/fnc/chapter$chap"
section_files = filter(endswith(".md"), readdir(dir, join=true))

for file in section_files[2:2]
    text = readlines(file, keep=true)
    start = findfirst(startswith("## Exercises"), text)
    isnothing(start) && continue
    mtch = findall(s -> !isnothing(match.(r"^[0-9]+\. ", s)), text[start+1:end]) .+ start
    push!(mtch, length(text) + 2)
    push!(text, "")
    push!(text, "")
    for i in 1:length(mtch)-1
        text[mtch[i]] = replace(text[mtch[i]], r"^[0-9]+\. " => "\n``````{exercise}\n")
        for j in mtch[i]+1:mtch[i+1]-1
            text[j] = lstrip(text[j], [' ', '\t'])
        end
        text[mtch[i+1]-1] = "``````\n" * text[mtch[i+1]-1]
    end
    for j in start+1:length(text)
        mtch = match(r"\s*\(([a-zA-Z\-]+)\)=", text[j])
        isnothing(mtch) && continue
        println("Found label in $file: $j, $(mtch.captures[1])")
        text[j] = "\n"
        jj = findnext(startswith("\n``````{exercise}"), text, j)
        text[jj] = replace(text[jj], "{exercise}\n" => "{exercise}\n:label: " * mtch.captures[1] * "\n")
    end

    open(file, "w") do file
        write.(file, text)
    end
end
