(chapter-localapprox)=

# 5. Piecewise interpolation

```{index} Yoda, The Empire Strikes Back
```

```{epigraph}
You must feel the Force around you. Here, between you...me...the tree...the rock...everywhere!

— Yoda, *The Empire Strikes Back*
```

In many scientific problems the solution is a function. Accordingly, our next task is to represent functions numerically. This task is more difficult and complicated than the one we faced in representing real numbers. With numbers, it's intuitively clear how one real value can stand for a small interval around it. But designating representatives for sets of functions is less straightforward—in fact, it's one of the core topics in computing. The process of converting functions into numerical representations of finite length is known as *discretization*.

Once we have selected a method of discretization, we can define numerical analogs of our two favorite operations on functions, differentiation and integration. These are linear operations, so the most natural numerical analogs are linear operations too. As we will see in many of the chapters following this one, a lot of numerical computing boils down to converting calculus to algebra, with discretization as the link between them.
