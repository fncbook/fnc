chap = 4
dir = "/Users/driscoll/Documents/GitHub/fnc/chapter$chap"
code_files = [
    "/Users/driscoll/Documents/GitHub/fnc/julia/chapter$chap.md",
    "/Users/driscoll/Documents/GitHub/fnc/matlab/chapter$chap.md",
    "/Users/driscoll/Documents/GitHub/fnc/python/chapter$chap.md"
]

repls = [
"Section $chap.1" => "$chap.1 @section-nonlineqn-rootproblem",
"Section $chap.2" => "$chap.2 @section-nonlineqn-fixed-point",
"Section $chap.3" => "$chap.3 @section-nonlineqn-newton",
"Section $chap.4" => "$chap.4 @section-nonlineqn-secant",
"Section $chap.5" => "$chap.5 @section-nonlineqn-newtonsys",
"Section $chap.6" => "$chap.6 @section-nonlineqn-quasinewton",
"Section $chap.7" => "$chap.7 @section-nonlineqn-nlsq",
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
