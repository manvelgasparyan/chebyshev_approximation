function [coefficients, polynomial, x_values, y_values] = Chebyshev_approximation (Input, N, m, a, b)
%==========================================================================
   %check whether the main input is a function handle or a dataset
   if isa (Input, 'function_handle') %If the input is a function handle
      f = Input;
      [coefficients, polynomial, x_values, y_values] = Chebyshev_function_approximation (f, N, m, a, b);
   else %If the input is a dataset
      Z = Input;
      [coefficients, polynomial, x_values, y_values]  = Chebyshev_data_approximation (Z, N, m);
   end
%==========================================================================
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coefficients, polynomial, x_values, y_values] = Chebyshev_function_approximation (f, N, m, a, b)
%==========================================================================
%==========================================================================
%==========================================================================
   %The first N Chebyshev coefficients of f over the interval [a, b];
   alpha = Chebyshev_coefficients_f (f, N, a, b);
   coefficients = alpha;
   %-----------   
   %The points of interest on the interval [a,b]
   interval = a:(b-a)/m:b;
   N_interval = numel(interval);
   x_values = interval;
   %----------
   y_values = zeros(N_interval, 1);
   %----------
   for i = 1:N_interval
    interval_i = interval(i);
    for n = 0:N
        c_n = @(x)Chebyshev_polynomial (x, n, a, b);
        c_in = c_n(interval_i);
        alpha_n = alpha(n+1);
        y_values(i) = y_values(i) + alpha_n*c_in;
    end
   end
   %----------   
   %Determine the explicit form of the polynomial
   syms x
   p_N = @(t)Chebyshev_series (t, coefficients, N, a, b);
   polynomial = p_N(x);
   polynomial = simplify(polynomial);   
%==========================================================================
%==========================================================================
%==========================================================================
%Computing absolute values of the the coefficients
coefficients = abs(coefficients);
%-----------
%Plotting the coefficients
figure(1)
subplot(1, 2, 1)
plot(0:1:N, coefficients, 'x', 'linewidth', 1.5, 'color', 'r')
   %------
   title_label = 'Coefficients of the Chebyshev polynomial of $f$';
   x_label = 'Indexing';
   y_label = 'Absolute values of the coefficients';
   my_legend = 'Chebyshev coefficients';
   %------
   set(gca, 'color', [0.75 0.75 0.75]) 
   grid on
   ax = gca;
   ax.GridColor = 'white';
   ax.GridAlpha = 0.2;
   grid minor
   ax.MinorGridColor = 'white';
   legend(my_legend, 'Location', 'Best');
   xlabel(x_label, 'fontweight','bold', 'Fontsize', 15, 'interpreter','latex')
   ylabel(y_label, 'fontweight','bold', 'Fontsize', 15, 'interpreter', 'latex')
   title(title_label, 'fontweight','bold', 'Fontsize', 15, 'interpreter', 'latex')
   lhs = min(coefficients) - 0.05*min(coefficients);
   rhs = max(coefficients) - 0.05*max(coefficients);
   ylim([lhs rhs])
%-----------
%The values of the original function f in the interval [a, b] using m+1
%division points
F = f(x_values);
%-----------
%plotting the approcimation
figure(1)
subplot(1,2,2)
%-----------
plot(x_values, F, 'linewidth', 1.5, 'color', 'r')
hold on
fplot(polynomial, [a b], '--', 'linewidth', 1.5, 'color', 'b')
%plot(x_values, y_values, '--', 'linewidth', 1.5, 'color', 'b')
hold off
   %------
   title_label = 'Approximation of $f$ using the first $N$ Chebyshev polynomials';
   x_label = 'The interval of interest';
   y_label = 'The original function and the Chebyshev approximation';
   my_legend = {'Original function', 'Chebyshev approximation'};
   %------
   set(gca, 'color', [0.75 0.75 0.75]) 
   grid on
   ax = gca;
   ax.GridColor = 'white';
   ax.GridAlpha = 0.2;
   grid minor
   ax.MinorGridColor = 'white';
   legend(my_legend, 'Location', 'Best');
   xlabel(x_label, 'fontweight','bold', 'Fontsize', 15, 'interpreter','latex')
   ylabel(y_label, 'fontweight','bold', 'Fontsize', 15, 'interpreter', 'latex')
   title(title_label, 'fontweight','bold', 'Fontsize', 15, 'interpreter', 'latex')
%==========================================================================
%==========================================================================
%==========================================================================
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coefficients, polynomial, x_values, y_values]  = Chebyshev_data_approximation (Z, N, m)
%----------
X = Z(:,1);
Y = Z(:, 2);
%----------
%The best fitting Chebyshev coefficients
alpha_hat = Chebyshev_coefficients_d (X, Y, N);
coefficients = abs(alpha_hat);
%----------
%The interval of intterest
a = X(1);
b = X(end);
%----------
%Iteration of the interval
x_values = a:(b-a)/m:b;
N_interval = numel(x_values);
%----------
%The Chebyshev series
p_N = @(t)Chebyshev_series (t, alpha_hat, N, a, b);
%----------
%The output is the predicted values
y_values = zeros(N_interval, 1);
for i = 1:N_interval
    x_values_i = x_values(i);
    y_values(i) = p_N(x_values_i);
end
%-----------
%Determine the explicit form of the polynomial
syms x
polynomial = p_N(x);
polynomial = simplify(polynomial);
%-----------
%Plotting the coefficients
figure(2)
%-----------
subplot(1, 2, 1)
plot(1:1:N+1, coefficients, 'x', 'linewidth', 1.5,  'color', 'k')
%-----------
   title_label = 'The absolute values of the best fitting Chebyshev coefficients';
   x_label = 'Indexing';
   y_label = 'Chebyshev coefficients';
   my_legend = 'Coefficients';
   %------
   set(gca, 'color', [0.75 0.75 0.75]) 
   grid on
   ax = gca;
   ax.GridColor = 'white';
   ax.GridAlpha = 0.2;
   grid minor
   ax.MinorGridColor = 'white';
   legend(my_legend, 'Location', 'Best');
   xlabel(x_label, 'fontweight','bold', 'Fontsize', 15, 'interpreter','latex')
   ylabel(y_label, 'fontweight','bold', 'Fontsize', 15, 'interpreter', 'latex')
   title(title_label, 'fontweight','bold', 'Fontsize', 15, 'interpreter', 'latex')
%-----------
subplot(1,2,2)
%plot(x_values, y_values, 'linewidth', 1.5, 'color', 'b')
fplot(polynomial, [a b], 'linewidth', 1.5, 'color', 'b')
hold on
plot(X, Y, 'x', 'linewidth', 1.5, 'color', 'r')
hold off
%-----------
   title_label = 'Data fitting using the Chebyshev polynomials';
   x_label = 'X Values';
   y_label = 'Comparison of the provided and predicted Y values';
   my_legend_2 = 'Available data';
   my_legend_1 = 'Predicted values';
   %------
   set(gca, 'color', [0.75 0.75 0.75]) 
   grid on
   ax = gca;
   ax.GridColor = 'white';
   ax.GridAlpha = 0.2;
   grid minor
   ax.MinorGridColor = 'white';
   legend(my_legend_1, my_legend_2, 'Location', 'Best');
   xlabel(x_label, 'fontweight','bold', 'Fontsize', 15, 'interpreter','latex')
   ylabel(y_label, 'fontweight','bold', 'Fontsize', 15, 'interpreter', 'latex')
   title(title_label, 'fontweight','bold', 'Fontsize', 15, 'interpreter', 'latex')
%-----------
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output] = Chebyshev_coefficients_f (f, N, a, b)
%==========
output = zeros(N+1, 1);
%==========
%The transformation
tau = @(t)(2*(t-a)./(b-a) - 1);
%==========
%The transformed weight function
u = @(t)(1./sqrt(1 - tau(t).^2));
%==========
%The first coefficient
C_0 = 1;
f_0 = @(t)f(t).*C_0.*u(t);
output(1) = (2/((b-a)*pi))*integral(f_0, a, b);
%==========
%The nth coefficient for n = 1,...,N
for n = 1:N
    C_n = @(t)Chebyshev_polynomial (t, n, a, b);
    f_n = @(t)f(t).*C_n(t).*u(t);
    output(n+1) = (4/((b - a)*pi))*integral(f_n, a, b);
end
%==========
output = output';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alpha_hat] = Chebyshev_coefficients_d (X, Y, N)
%==========
%The objective function of minimization
error = @(alpha)Chebyshev_error (alpha, X, Y, N);
%==========
%The provided guess points
alpha_guess = rand(N+1, 1);
%==========
%Determination of the best fitting Chebyshev coefficients
alpha_hat = fmincon(error, alpha_guess);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C_n] = Chebyshev_polynomial (x, n, a, b)
%==========
%The nth regular Chebyshev polynomial
T_n = @(t)chebyshevT (n, t);
%==========
%Transformation from [a,b] to [-1, 1]
x = 2*(x-a)./(b-a) - 1;
%==========
%The trnsformed nth Chebyshev polynomial
C_n = T_n(x);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p_N] = Chebyshev_series (t, alpha, N, a, b)
%==========
%The output is the Chebyshev series of degree N
p_N = 0;
%==========
%The nth coefficient for n = 0, 1,...,N
for n = 0:N
    C_n = Chebyshev_polynomial (t, n, a, b);
    p_N = p_N + alpha(n+1)*C_n;
end
%==========
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [error] = Chebyshev_error (alpha, X, Y, N)
%==========
%The length of the dataset
M = numel(X);
%==========
%The interval of interest
a = X(1);
b = X(M);
%==========
%The Chebyshev polynomial of degree N
p_N = @(t)Chebyshev_series (t, alpha, N, a, b);
%==========
%Predicted values using Chebyshev series
Chebyshev_predicted = zeros(M, 1);
for i = 1:M
    X_i = X(i);
    Chebyshev_predicted(i) = p_N(X_i); 
end
%==========
%The difference between available data and the corresponding predicted
%values using Chebyshev series
z = Chebyshev_predicted - Y;
%==========
%The error
error = z'*z;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
