---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
exec(open("../../../python/FNC_init.py").read())
```
[**Demo %s**](#demo-insight-image)

```{index} ! Julia; Images
```

We will use a test image from the well-known `scikit-image` package.

```{code-cell}
from skimage import data as testimages
img = getattr(testimages, "coffee")()
imshow(img)
```

The variable `img` is a matrix.

```{code-cell}
size(img)
```

However, its entries are colors, not numbers.

```{code-cell}
print(f"image has shape {img.shape}")
print(f"first pixel has value {img[0, 0]}")
```

The three values at each pixel are for intensities of red, green, and blue. We can convert each of those layers into an ordinary matrix of values between 0 and 255, which is maximum intensity.

```{code-cell}
R = img[:, :, 0]
print("upper left corner of the red plane is:")
print(R[:5, :5])
print(f"red channel values range from {R.min()} to {R.max()}")
```

It may also be convenient to convert the image to grayscale, which has just one layer of values from zero (black) to one (white).

```{code-cell}
from skimage.color import rgb2gray
A = rgb2gray(img)
A[:5, :5]
print("upper left corner of grayscale:")
print(A[:5, :5])
print(f"gray values range from {A.min()} to {A.max()}")
```

```{code-cell}
imshow(A, cmap='gray')
axis('off');
```

Some changes we make to the grayscale matrix are easy to interpret visually.

```{code-cell}
imshow(flipud(A), cmap='gray')
axis('off');
```
