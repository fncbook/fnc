using JSON

"""
    find_objects(data)

# Arguments
- `data`: The JSON data structure

# Returns
- `Vector{Dict}`: Array of dictionaries that have key=value
"""
function find_objects(data, key="output_type", value="execute_result")
    results = Vector{typeof(data)}()

    # Helper function to recursively traverse the data structure
    function traverse(obj)
        if isa(obj, JSON.Object)
            # Check if this dictionary has type="output"
            if haskey(obj, key) && obj[key] == value
                push!(results, obj)
            end

            # Recursively traverse all values in the dictionary
            for (key, value) in obj
                traverse(value)
            end

        elseif isa(obj, Array)
            # Recursively traverse all elements in the array
            for item in obj
                traverse(item)
            end
        end
        # For primitive types (String, Number, Bool, etc.), do nothing
    end

    traverse(data)
    return results
end

# Validate that HTML build files have output content.
# root = "separate/python"
root = "."
root = joinpath(root, "_build", "html")
# last_in_chapter = ["stability", "structure", "house", "nlsq", "adaptive", "zerostability", "dimreduce", "precond", "improper", "galerkin", "boundaries", "wave", "nonlinear-1"]
found = Dict()
for fname in filter(endswith(".json"), readdir(root))
    println("Checking $fname...")
    json = JSON.parsefile(joinpath(root, fname))
    for lang in ["julia", "python", "matlab"]
        for obj in find_objects(json, "sync", lang)
            for er in find_objects(obj, "output_type", "error")  # for python and julia
                get!(found, fname, Dict())[lang] = er["evalue"]
            end
            for er in find_objects(obj, "name", "stderr")  # for matlab
                if startswith(get!(er, "text", "") , "Error")
                    get!(found, fname, Dict())[lang] = er["text"]
                end
            end
        end
    end
end

intended = Dict("matrices.json" => 3, "structure.json" => 3, "linear-systems.json" => 2)

println("\n\nREPORT\n------")
for (key, value) in found
    println("File: $key")
    if haskey(intended, key)
        if length(value) == intended[key]
            println("  Found expected number of errors: $(length(value))")
            continue
        end
    end
    for (lang, err) in value
        println("  $lang error: $err")
    end
    println()
end
