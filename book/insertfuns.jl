using FundamentalsNumericalComputation

fncdir = dirname(pathof(FNC))

for root in ["intro","linsys","leastsq","nonlineqn","localapprox","ivp","appendix"]
	for mdfile in filter(s->endswith(s,".md"),readdir(root))
		str = readlines(joinpath(root,mdfile))
		itext = findfirst(occursin.("{proof:function}",str))
		offset = 0
		while !isnothing(itext)
			itext += offset
			funcname = match(r"`+{proof:function} (.*)",str[itext]).captures[1]
			@show root,mdfile,itext,funcname

			funcstr = ""

			for srcfile in filter(s->startswith(s,"chapter"),readdir(fncdir))
				src = readlines(joinpath(fncdir,srcfile))
				isrc = findfirst(occursin.("function $funcname",src))
				isnothing(isrc) && continue

				k = isrc + findfirst(occursin.("return",src[isrc+1:end]))
				endsrc = k + findfirst(startswith.(src[k+1:end],"end"))

				k = findlast(startswith.(src[1:isrc-1],"\"\"\""))
				startsrc = findlast(startswith.(src[1:k-1],"\"\"\""))

				funcstr = src[startsrc:endsrc]
				break
			end

			n = itext + findfirst(occursin.("{code-block}",str[itext+1:end]))
			offset = n + findfirst(startswith.(str[n+1:end],"```"))

			newstr = [":lineno-start: 1";funcstr]
			str = [str[1:n];newstr;str[offset:end]]
			offset = n + length(newstr) + 1
			
			itext = findfirst(occursin.("{proof:function}",str[offset+1:end]))
		end

		open(joinpath(root,"$mdfile"),"w") do io 
			[ println(io,s) for s in str ];
		end
	end
end