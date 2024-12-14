(section-leastsq-next)=
# Next steps

The least-squares problem has been widely studied and used, and only seems to become more important in this era of ever-increasing amounts of data.  A good reference for numerical methods is the monograph by Björck {cite}`bjorckNumericalMethods1996`.  Some theoretical results can be found in Higham {cite}`highamAccuracyStability2002`; a brief and advanced discussion can be found in Golub and Van Loan {cite}`golubMatrixComputations1996`.

Note that a vast literature can also be found in statistics for what is referred to as **data regression**, or simply **regression**. Nonlinear methods for least-squares fitting of data will be discussed in the following chapter.

In modern applications one may have to deal with so-called online fitting, in which new data must continually be incorporated with old. More recent sources address related issues, e.g., in  Hansen, Pereyra, and Scherer {cite}`hansenLeastSquares2013` and in Teunissen {cite}`teunissenDynamicData2001`.  The problem of geodesy and GPS positioning are discussed in some detail in Strang and Borre {cite}`strangLinearAlgebra1997`; for these applications, they describe how the updating of least squares leads to Kalman filtering.
