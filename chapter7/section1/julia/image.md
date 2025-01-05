---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
```{code-cell}
include("../../../julia/FNC_init.jl")
```
[**Demo %s**](#demo-insight-image)

```{index} ! Julia; Images
```

The `Images` package has many functions for image manipulation, and `TestImages` has some standard images to play with.

```{code-cell}
using Images, TestImages
img = testimage("mandrill")
```

The variable `img` is a matrix.

```{code-cell}
size(img)
```

However, its entries are colors, not numbers.

```{code-cell}
img[100, 10]
```

```{index} ! Julia; eltype
```

You can use `eltype` to find out the type of the elements of any array.

```{code-cell}
eltype(img)
```

It's possible to extract matrices of red, green, and blue intensities, scaled from 0 to 1.

```{code-cell}
R, G, B = red.(img), green.(img), blue.(img);
@show minB, maxB = extrema(B);
```

Or we can convert the pixels to gray, each pixel again scaled from 0 to 1.

```{code-cell}
Gray.(img)
```

In order to do our usual operations, we need to tell Julia that we want to interpret the elements of the image matrix as floating-point values.

```{code-cell}
A = Float64.(Gray.(img))
A[1:4, 1:5]
```

We can use `Gray` to reinterpret a matrix of floating-point values as grayscale pixels.

```{code-cell}
Gray.(reverse(A, dims=1))
```
