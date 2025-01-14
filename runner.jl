using JSON

# chapter = 1
# section = 1
lang = "julia"

function editor(cell)
    if cell["cell_type"] == "code" && haskey(cell["metadata"], "tags")
        tags = cell["metadata"]["tags"]
        if "remove-cell" in tags
            return nothing
        elseif "remove-output" in tags
            cell["outputs"] = []
        elseif "remove-input" in tags
            cell["source"] = []
        end
    end
    return cell
end


for chapter in 3:13, section in 1:9
    dir = joinpath(pwd(), "chapter$chapter/section$section/$lang/")
    nbs = try
        filter(x -> endswith(x, "nbconvert.ipynb") , readdir(dir))
    catch
        continue
    end
    for nb in nbs
        # cmd = `/Users/driscoll/mambaforge/envs/myst/bin/jupyter-nbconvert --to notebook --execute $file`
        # run(cmd)

        contents = open(joinpath(dir, nb)) do io
            JSON.parse(io)
        end

        # Julia ONLY
        if lang == "julia"
            contents["cells"][1]["metadata"]["tags"] = ["remove-cell"]
        end
        original = replace(nb, ".nbconvert.ipynb" => ".ipynb")
        mkpath(joinpath(dir, "originalnb"))
        open(joinpath(dir, "originalnb", original), "w") do io
            JSON.print(io, contents, 4)
        end

        contents["cells"] = filter(!isnothing, editor.(contents["cells"]))
        open(joinpath(dir, original), "w") do io
            JSON.print(io, contents, 4)
        end

    end
end


##

for (root, dirs, files) in walkdir("chapter13"), f in files
    if contains(root, "originalnb") && endswith(f, ".nbconvert.ipynb")
        newf = replace(f, ".nbconvert.ipynb" => ".ipynb")
        mv(joinpath(root, f), joinpath(root, newf))
        # println(joinpath(root, f), " to  ", joinpath(root, newf))
    end
end
