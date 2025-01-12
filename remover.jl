using JSON
file = "chapter3/section1/matlab/pirate.ipynb"
file = "chapter1/section1/matlab/accuracy.ipynb"
function isdoomed(cell)
    return haskey(cell, "metadata") && haskey(cell["metadata"], "tags") && "remove-cell" in cell["metadata"]["tags"]
end
j = open(file) do io
    JSON.parse(io)
end

deleteat!(j["cells"], isdoomed.(j["cells"]))

open("tmp.ipynb", "w") do io
    JSON.print(io, j)
end
