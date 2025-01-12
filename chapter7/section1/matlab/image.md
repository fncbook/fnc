---
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
cd  /Users/driscoll/Documents/GitHub/fnc/matlab
FNC_init
```
[**Demo %s**](#demo-insight-image)

```{index} ! Julia; Images
```

MATLAB ships with a few test images to play with.

```{code-cell}
A = imread('peppers.png');
color_size = size(A)
```

Use `imshow` to display the image.

```{code-cell}
imshow(A)
```

The image has three layers or channels for red, green, and blue. We can deal with each layer as a matrix, or (as below) convert it to a single matrix indicating shades of gray from black (0) to white (255). Either way, we have to explicitly convert the entries to floating-point values rather than integers. 


```{code-cell}
A = im2gray(A);   % collapse from 3 dimensions to 2
gray_size = size(A)
imshow(A)
```

Before we can do any numerical computation, we need to convert the image to a matrix of floating-point numbers.

```{code-cell}
A = double(A);
```
