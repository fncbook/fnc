chap = 1

dir = "/Users/driscoll/Documents/GitHub/fnc/chapter$chap"
section_files = filter(endswith(".md"), readdir(dir, join=true))

for file in section_files[1:end]
    text = readlines(file, keep=true)
    i = 0
    while i < length(text)
        i += 1
        mtch = match(r"^(:+){math}", text[i])
        isnothing(mtch) && continue
        println("Found math fence in $file: $i, $(mtch.captures[1])")
        n = length(mtch.captures[1])
        text[i] = replace(text[i], ":" => "`")
        i += 1
        while i <= length(text) && !startswith(text[i], ":"^n)
            i += 1
        end
        text[i] = replace(text[i], ":" => "`")
    end
    open(file, "w") do file
        write.(file, text)
    end
end
