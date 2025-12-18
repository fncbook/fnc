using Printf

function get_code_blocks(chap, lang)
    header =
    if lang == "python"
        yaml = """
kernelspec:
  display_name: Python 3
  language: python
  name: python3
"""
        init = readlines("/Users/driscoll/Documents/GitHub/fnc/python/FNC_init.py")
        (; yaml, init)
    elseif lang == "matlab"
        yaml = """
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
"""
        init = replace.(readlines("/Users/driscoll/Documents/GitHub/fnc/matlab/FNC_init.m"), "FNC-matlab/" => "../FNC_matlab/")
        (; yaml, init)
    elseif lang == "julia"
        yaml = """
        kernelspec:
          display_name: Julia 1
          language: julia
          name: julia-1.12
        """
        init = readlines("/Users/driscoll/Documents/GitHub/fnc/julia/FNC_init.jl")
        (; yaml, init)
    end

    # read functions and demos contents
    code_file = readlines("/Users/driscoll/Documents/GitHub/fnc/$lang/chapter$chap.md")
    func_start = something(findfirst(contains("## Functions"), code_file), length(code_file) + 1)
    example_start = something(findfirst(contains("## Examples"), code_file), length(code_file) + 1)

    blocks = Dict()
    idx = func_start + 1
    while idx < example_start
        fun_idx = findnext(contains(r"\(function-.*\)"), code_file, idx)
        @show fun_idx
        isnothing(fun_idx) && break
        tag = match(r"\((.*-.*)\)", code_file[fun_idx]).captures[1]
        drop_idx = findnext(contains("{dropdown}"), code_file, fun_idx)
        backticks = match(r"(`+){dropdown}", code_file[drop_idx]).captures[1]
        drop_end = findnext(contains(Regex("$backticks")), code_file, drop_idx+1)
        # blocks[tag] = code_file[drop_idx:drop_end]
        blocks[tag] = code_file[drop_idx+2:drop_end-1]
        idx = drop_end + 1
    end

    idx = example_start + 1
    while idx < length(code_file)
        exam_idx = findnext(contains(r"(demo-.*-)"), code_file, idx)
        isnothing(exam_idx) && break
        tag = match(r"\((.*)\)", code_file[exam_idx]).captures[1]
        println("Found demo tag: $tag")
        drop_idx = findnext(contains("{dropdown}"), code_file, exam_idx)
        backticks = match(r"(`+){dropdown}", code_file[drop_idx]).captures[1]
        drop_end = findnext(contains(Regex("$backticks")), code_file, drop_idx+1)
        # blocks[tag] = code_file[drop_idx:drop_end]
        blocks[tag] = code_file[drop_idx+2:drop_end-1]
        idx = drop_end + 1
    end
    return header, blocks
end

function transfer_content(dir, new_dir, chap, lang)
    header =  (; yaml=[], init=[])
    blocks = Dict()
    if chap > 0
        println("Getting code blocks for chapter $chap")
        header, blocks = get_code_blocks(chap, lang)
    end

    # fix up code includes
    if lang == "python"
        try
            fn = @sprintf("chapter%02d.py", chap)
            cp(joinpath("/Users/driscoll/Documents/GitHub/fnc/python/fncbook/fncbook", fn), joinpath(new_dir, fn); force=true)
        catch
        end
        for (key, val) in pairs(blocks)
            blocks[key] = replace.(val, "fncbook/fncbook/" => "")
        end
    elseif lang == "julia"
        try
            fn = @sprintf("chapter%02d.jl", chap)
            cp(joinpath("/Users/driscoll/Documents/GitHub/fnc/julia/FNCFunctions/src", fn), joinpath(new_dir, fn); force=true)
        catch
        end
        for (key, val) in pairs(blocks)
            blocks[key] = replace.(val, "FNCFunctions/src/" => "")
        end
    elseif lang == "matlab"
        cp("/Users/driscoll/Documents/GitHub/fnc/matlab/FNC-matlab", "/Users/driscoll/Documents/GitHub/fnc/separate/matlab/FNC_matlab"; force=true)
        for (key, val) in pairs(blocks)
            blocks[key] = replace.(val, "FNC-matlab/" => "../FNC_matlab/")
        end
    end

    # transfer other files
    files = filter(!endswith(".md"), readdir(dir))
    foreach(file -> cp(joinpath(dir, file), joinpath(new_dir, file); force=true), files)

    # process markdown files
    files = filter(endswith(".md"), readdir(dir))
    excerpt = []
    for file in files
        println("\nProcessing $file")
        file_path = joinpath(dir, file)
        md_content = readlines(file_path)
        open(joinpath(new_dir, file), "w") do f
            # write out the header
            idx = findnext(startswith("---"), md_content, 2)
            if !isnothing(idx) && chap > 0
                write(f, "---\n")
                foreach(line -> write(f, line * "\n"), md_content[2:idx-1])
                write(f, header.yaml * "---\n")
                write(f, "```{code-cell}\n:tags: [remove-cell]\n")
                foreach(line -> write(f, line * "\n"), header.init)
                write(f, "```\n\n")
                idx += 1
            else
                idx = 1
            end

            # work up to each tab set
            while idx < length(md_content)
                tab_idx = findnext(contains(r"(`+){tab-set}"), md_content, idx)
                isnothing(tab_idx) && break
                foreach(line -> write(f, line * "\n"), md_content[idx:tab_idx-1])
                backticks = match(r"(`+){tab-set}", md_content[tab_idx]).captures[1]
                tab_end = findnext(contains(Regex("$backticks")), md_content, tab_idx+1)
                excerpt = md_content[tab_idx:tab_end]
                @show tab_idx, tab_end

                # find the tab-item for the desired language
                item_idx = findfirst(contains(Regex("(`+){tab-item}.*$lang")), lowercase.(excerpt))
                item_ticks = match(r"(`+){tab-item}", excerpt[item_idx]).captures[1]
                item_end = findnext(contains(Regex("$item_ticks")), excerpt, item_idx + 1)
                # sync_idx = findfirst(contains(":sync: $lang"), lowercase.(excerpt))
                # isnothing(sync_idx) && error("File $file does not contain :sync: at $tab_idx")
                @show item_idx, item_end
                embed_idx = findnext(contains(r"{embed}"), excerpt, item_idx + 1)
                if isnothing(embed_idx)
                    println("no embed found")
                    # just include raw content
                    foreach(excerpt[item_idx+1:item_end-1]) do line
                        !startswith(line, ":sync:") && write(f, line * "\n")
                    end
                else
                    println("embed at $embed_idx")
                    embed = match(r".*{embed}[ ]*#(.*)", excerpt[embed_idx])
                    tag = embed.captures[1]
                    foreach(line -> write(f, line * "\n"), blocks[tag])
                end
                idx = tab_end + 1
            end

            # write out the rest of the file
            foreach(line -> write(f, line * "\n"), md_content[idx:end])
        end
    end
end

##

lang = "matlab"
for chap in 1:13
    println("\nProcessing chapter $chap for $lang")
    dir = "/Users/driscoll/Documents/GitHub/fnc/chapter$chap"
    new_dir = "/Users/driscoll/Documents/GitHub/fnc/separate/$lang/chapter$chap"
    mkpath(new_dir)
    mkpath(new_dir * "/figures")
    transfer_content(dir, new_dir, chap, lang)
end

println("Processing appendix for $lang")
dir = "/Users/driscoll/Documents/GitHub/fnc/appendix"
new_dir = "/Users/driscoll/Documents/GitHub/fnc/separate/$lang/appendix"
mkpath(new_dir)
transfer_content(dir, new_dir, 0, lang)

##

dir = "/Users/driscoll/Documents/GitHub/fnc"
new_dir = "/Users/driscoll/Documents/GitHub/fnc/separate/$lang"
for fn in ["home.md", "genindex.md", "refs.md", "FNC.bib", "_static", "frontmatter"]
    cp(joinpath(dir, fn), joinpath(new_dir, fn); force=true)
end
cp(joinpath(dir, "$lang/setup.md"), joinpath(new_dir, "setup.md"); force=true)

if lang == "python"
    cp("/Users/driscoll/Documents/GitHub/fnc/python/roswelladj.mat", "/Users/driscoll/Documents/GitHub/fnc/separate/$lang/chapter8/roswelladj.mat"; force=true)
    cp("/Users/driscoll/Documents/GitHub/fnc/python/voting.mat", "/Users/driscoll/Documents/GitHub/fnc/separate/$lang/chapter7/voting.mat"; force=true)
elseif lang == "julia"
    cp("/Users/driscoll/Documents/GitHub/fnc/julia/roswell.jld2", "/Users/driscoll/Documents/GitHub/fnc/separate/$lang/chapter8/roswell.jld2"; force=true)
    cp("/Users/driscoll/Documents/GitHub/fnc/julia/smallworld.jld2", "/Users/driscoll/Documents/GitHub/fnc/separate/$lang/chapter8/smallworld.jld2"; force=true)
    cp("/Users/driscoll/Documents/GitHub/fnc/julia/voting.jld2", "/Users/driscoll/Documents/GitHub/fnc/separate/$lang/chapter7/voting.jld2"; force=true)
elseif lang == "matlab"
    for chap in [4, 6, 12, 13]
        for fn in filter(contains(Regex("f$chap.*\\.m")), readdir("/Users/driscoll/Documents/GitHub/fnc/matlab"))
            cp(joinpath("/Users/driscoll/Documents/GitHub/fnc/matlab", fn), joinpath("/Users/driscoll/Documents/GitHub/fnc/separate/matlab/chapter$chap", fn); force=true)
        end
    end
    for chap in 2:13
        cp("/Users/driscoll/Documents/GitHub/fnc/matlab/redsblues.m", "/Users/driscoll/Documents/GitHub/fnc/separate/matlab/chapter$chap/redsblues.m"; force=true)
    end
    cp("/Users/driscoll/Documents/GitHub/fnc/matlab/roswelladj.mat", "/Users/driscoll/Documents/GitHub/fnc/separate/$lang/chapter8/roswelladj.mat"; force=true)
    cp("/Users/driscoll/Documents/GitHub/fnc/matlab/smallworld.mat", "/Users/driscoll/Documents/GitHub/fnc/separate/$lang/chapter8/smallworld.mat"; force=true)
    cp("/Users/driscoll/Documents/GitHub/fnc/matlab/voting.mat", "/Users/driscoll/Documents/GitHub/fnc/separate/$lang/chapter7/voting.mat"; force=true)
end
