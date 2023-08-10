# Chebyshev Function and Data Approximation - MATLAB


This GitHub repository contains two MATLAB functions for approximating functions and datasets using Chebyshev polynomials, respectively. Chebyshev polynomials are a set of orthogonal polynomials commonly used for numerical approximation and interpolation. For a deeper understanding of Chebyshev polynomials, please refer to [1].

The first function, named "Chebyshev_function_approximation", can be used to approximate a continuous function f on any interval '[a, b]' with the partial sum 'p_N' of the Chebyshev expansion of degree 'N'. It takes the following inputs:

'f' - The continuous function to be approximated. It should be defined as a function handle, such as 'f = @(x)-(3./2).*exp(-4.*x)'.
'a' - The lower bound of the interval.
'b' - The upper bound of the interval.
'N' - The degree of the partial sums of the Chebyshev expansion.
'm' - The number of points on the interval '[a, b]'.
The function provides the following outputs:

polynomial - The polynomial 'p_N' of degree 'N' that approximates the function 'f' on the interval '[a, b]'.
coefficients - The coefficients of the Chebyshev functions in the expression of the polynomial.
x_values - The points on the interval '[a, b]' at which the polynomial 'p_N' is evaluated.
y_values - The values of the Chebyshev polynomial 'p_N' at 'm' points of the interval '[a, b]'.
Visualization of the absolute values of the Chebyshev coefficients.
Comparison plots between the original function 'f' and the polynomial 'p_N,' which approximates 'f' over the interval [a, b]."
Example:

f = @(x)-(3./2).*exp(-4.*x);

N = 10;

m = 100;

a = -2;

b = 4;

[coefficients, polynomial, x_values, y_values] = Chebyshev_function_approximation (f, N, a, b, m);
