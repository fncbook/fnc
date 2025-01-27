# Preface to the Julia edition

The invention of MATLAB introduced a new paradigm within research of numerical computation. Those concerned primarily with prototyping and perfecting algorithms, particularly those involving lots of linear algebra, optimization, and differential equations, were happy to adopt MATLAB as a primary computing environment. It offered concise syntax for such problems, freedom from variable types and declarations, compiling, and linking programs, convenient tools for analyzing results, and cross-platform uniformity. There were drawbacks, however, when it came to performance, scalability, and language features beyond the manipulation of vectors and matrices. While MATLAB has steadily made serious progress on closing the performance gap and introducing new language features, there remain computing tasks for which it is not ideally suited.

The landscape changed when the SciPy, NumPy, and Matplotlib packages for scientific computing in Python became stable and polished. These enabled Python to offer a fully featured MATLAB alternative that is free, open, and tightly integrated with a much larger world of computing. Among these valuable advances, though, a serious compromise lurked: performance got worse, often by orders of magnitude. While there are notable efforts to overcome the bottlenecks for many use cases, Python continues to face deep and steep performance challenges as a general-purpose scientific computing ecosystem.

Julia was designed from its inception to prioritize numerical scientific computing. It has reaped the benefits of learning from decisions and adaptations made in MATLAB and Python, borrowing the best parts from them and tackling their deficiencies. Julia's older cousins enjoy a big head start, so it's impossible to know what the size of Julia's niche will ultimately be, but interest has continued to build.

Why teach using Julia? The immediate benefits of Julia over MATLAB for the material in this text include:

* Julia allows Unicode characters, such as Greek letters, subscripts, and symbols, as variable names and operators, which makes code look more like mathematics.
* Julia makes it effortless to define functions inside of scripts as well as other functions.
* Julia's broadcasting syntax clarifies how to apply functions elementwise to arrays.
* Comprehensions are convenient and concise ways to construct vectors and matrices.
* Julia makes it easier to define keyword and optional function arguments.

There are also differences that cut both ways; for instance, Julia is often stricter about data types and sizes, which makes it more verbose and more prone to error, but arguably less likely to finish with unexpected results. 

Moreover, there are some tradeoffs. MATLAB ships with and installs everything needed for this text, while Julia requires a small installation effort to get started. MATLAB's documentation is superior, and it's easier to get accurate help on the Internet. MATLAB's integrated desktop, particularly the debugger, are not yet fully matched in Julia.
        
There is also a wider context to consider. Julia skills are more likely to be directly applicable in, or more easily transferrable to, high-performance applications. Julia interoperates easily with Python, R, C, and even MATLAB. Julia is native to the widely used Jupyter notebook system that currently dominates data science. Not least, as a free and open-source environment, Julia enables fully reproducible computing, which is increasingly appreciated as essential to long-term progress in research.

## What to expect from Julia

Unlike MATLAB and Python, Julia is just-in-time (JIT) compiled, not interpreted. As a result, large packages, including several supporting the code in this book, can take a few seconds to load.[^version] Furthermore, if you make a change to one of your own functions, or apply it to new types of function arguments, Julia may hesitate a moment while it compiles the necessary code. On slower hardware, or with frequent revision, the lag can become irritating.

[^version]: It is strongly recommended that you use at least version 1.6 of Julia. In earlier versions, the wait to compile and load packages can become long.

## What to expect from this book

*Supplemental material, including animations, downloadable code and examples, suggested projects, and more, can be found at https://bookstore.siam.org/ot177/bonus.

We do not attempt to teach how to become a great Julia programmer. That goal is too ambitious when stacked alongside the mathematical ones. Instead we hope to exhibit decent style and avoid promoting bad habits. But when there is a conflict, clarity and simplicity usually override performance concerns.[^fast]  A virtue of Julia is that one can start with a working straightforward code that can be adapted and improved to meet performance demands; we are introducing just the first stage of this process.

[^fast]: In a fast-changing language like Julia, yesterday's performance roadblocks can disappear anyway.

One choice advanced users might question is that we address vectors and matrices as starting from index 1, rather than using more general constructs such as `eachindex`, `begin`, and `first`. Our mathematics makes those specific references too, and in most languages, one must learn to deal directly with the difference between, say, 1-indexing and 0-indexing. 

Another notable choice we have made is the use of the popular `Plots` package for graphics. There are many other fine choices, including and not limited to `PyPlot`, `Makie`, and `PlotlyJS`, but we needed to be concrete.

Beyond that, we touch briefly on available packages that offer advanced functionality for the problem types we study. Our hope is that the student will not only learn fundamentals by working with simple codes in the book, but also learn about the existence and syntax of some power tools for serious applications.

Finally, we have made no mention of one of Julia's defining features, *multiple dispatch*. Using its power wisely edges into advanced software design and usually requires a bird's-eye view of problems that beginners lack. We want to let students expend as much of their cognitive budget as possible on the mathematical principles that have universal application. 

## Acknowledgments

We are indebted to Qinying Chen, Hugo Diaz, Mary Gockenbach, Aidan Hamilton, Pascal Kingsley Kataboh, Lindsey Jacobs, Ross Russell, and Jerome Troy, who made this text more accurate and more readable with their sharp eyes and great suggestions. And we are deeply grateful to Paula Callaghan at SIAM, whose patience, dedication, and wisdom were crucial to seeing this through to the end.
