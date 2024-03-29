Seeking to achieve greater flexibility in the use of test functions in the [Julia programming language](https://julialang.org/), we have initiated the development of a test function library in this language. This package, named `OptTestFun`, provides scalar functions defined in $\mathbb{R}^n$, along with their gradient vectors and Hessian matrices. Its aim is to study the performance of algorithms in the context of Continuous Optimization. We are including functions of various natures in the library, such as quadratic functions, functions involving trigonometric, exponential, and logarithmic expressions, as well as functions with sums of squares. Some of these functions were chosen from the website [`http://www.sfu.ca/~ssurjano`](http://www.sfu.ca/~ssurjano).

To install OptTestFun, you can access the Julia terminal and, in the package management mode, type the command add `https://github.com/AllanAmorim/OptTestFun.jl`. After installation, simply load the package using the command using OptTestFun. To access a specific objective function, along with its respective gradient vector and Hessian matrix, simply type the function name followed by `.obj` for the objective function, `.grad` for the gradient, and `.hess` for the Hessian in the terminal. The dimensionality of a function in OptTestFun is indirectly defined through the variable used. For example, the Rosenbrock function is implemented in this package and can be used via the command `rosen.obj(x)`, where `x` is an array of dimension $n$. Similarly, its gradient and Hessian matrix at a point $x$ can be used through the commands `rosen.grad(x)` and `rosen.hess(x)`, respectively.
