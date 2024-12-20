chap = 7
dir = "/Users/driscoll/Documents/GitHub/fnc/chapter$chap"
code_files = [
    "/Users/driscoll/Documents/GitHub/fnc/julia/chapter$chap.md",
    "/Users/driscoll/Documents/GitHub/fnc/matlab/chapter$chap.md",
    "/Users/driscoll/Documents/GitHub/fnc/python/chapter$chap.md"
]

repls = [
"Section $chap.1" => "$chap.1 @section-matrixanaly-insight",
"Section $chap.2" => "$chap.2 @section-matrixanaly-evd",
"Section $chap.3" => "$chap.3 @section-matrixanaly-svd",
"Section $chap.4" => "$chap.4 @section-matrixanaly-symm-eig",
"Section $chap.5" => "$chap.5 @section-matrixanaly-dimreduce",
]

for file in code_files
    code = readlines(file);
    for repl in repls
        code = replace.(code, repl)
    end
    open(file,"w") do file
        foreach(s->write(file, s*"\n"), code)
    end
end
