module OptTestFun

export rosen, dp, rast, sphere, trid, squares, stang

struct otf
    obj  :: Function
    grad :: Function
    hess :: Function 
end



# In this code, some of the expressions of the objective functions were adapted from
# the website https://www.sfu.ca/~ssurjano/optimization.html.



############### ROSENBROCK FUNCTION - DESCRIPTION ###############

# The Rosenbrock function, also known as the 'banana function', is a classic nonlinear optimization problem used for
# testing optimization algorithms. In $n$ dimensions, the Rosenbrock function is defined by the following expression:

# $$f(x)=\sum_{i=1}^{n-1} [100 (x_{i+1} - x_i^2)^2 + (x_i - 1)^2] $$ 

# This function consists of a weighted sum of quadratic and linear terms, where each term involves two adjacent variables. 
# The Rosenbrock function has the distinctive feature of having a global minimum at $(1, 1, ..., 1)$, 
# where the value of the function is zero. However, for $n > 2$, this function exhibits narrow and elongated valleys, 
# posing a challenge for many optimization algorithms. The function is often used to evaluate the effectiveness
# of optimization algorithms across a variety of applications.




############### ROSENBROCK FUNCTION (OBJECTIVE) ############### 

function rf(x::Vector{T}) where T <: Real # x is a vector of real coordinates
    n = length(x)   # Vector dimension
    sum = 0         # First part of the function
    for i in 1:length(x)-1 # For each element 'x[i]' of the vector 'x'...
        sum += 100 * (x[i+1] - x[i]^2)^2 + (1 - x[i])^2  # Part of the function involving the summation
    end
    return sum     # Sum of the first part with the second part
end


############### ROSENBROCK FUNCTION (GRADIENT) ############### 

function rg(x::Vector{T}) where T <: Real # x is a vector of real coordinates
    n = length(x)   # Vector dimension
    g = zeros(n)    # Creates a vector of n coordinates initialized to 0
    g[1] = 200 * (x[2] - x[1]^2) * (-2 * x[1]) + 2 * (x[1] - 1)   # Update the first coordinate of the vector with this result
    for i in 2:length(x)-1
        g[i] = 200 * (x[i] - x[i-1]^2) - 400 * x[i] * (x[i+1] - x[i]^2) + 2 * (x[i] - 1) # Update coordinates from 2 to n-1 with these results 
    end
    g[end] = 200.0 * (x[end] - x[end-1]^2)  # Update the last coordinate with this result
    return g
end

############### ROSENBROCK FUNCTION (HESSIAN) ############### 

function rh(x::Vector{T}) where T <: Real # x is a vector of real coordinates
    n = length(x)   # Vector dimension
    H = zeros(length(x), length(x))   # Creates a matrix of zeros with dimensions n x n

    # For the first row of the matrix:
    H[1, 1] = -400 * (x[2] - x[1]^2) + 800 * x[1]^2 + 2  # Replace the first result in the line with this value.
    H[1, 2] = -400.0 * x[1] # Replace the second result of the row with this value.

    # For rows 2 to n-1
    for i in 2:length(x)-1 
        H[i, i-1] = -400 * x[i-1]   # Replace the result of row i column i-1 with this value.
        H[i, i] = 200 - 400 * (x[i+1] - x[i]^2) + 800 * x[i]^2 + 2  # Replace the result of row i column i with this value.
        H[i, i+1] = -400.0 * x[i]  # Replace the result of row i column i+1 with this value.
    end

    # For the last row.
    H[end, end-1] = -400.0 * x[end-1]  # Replace the result of row n column n-1 with this value.
    H[end, end] = 200.0  # Replace the result of row n column n with this value.

    return H  
end






############### DIXON-PRICE FUNCTION - DESCRIPTION ###############

# The Dixon-Price function is a widely used optimization benchmark in testing optimization algorithms. 
# In $n$ dimensions, the Dixon-Price function is expressed as:

# $$ f(x) = (x_1 - 1)^2 + \sum_{i=2}^{n} i \cdot (2x_i^2 - x_{i-1})^2 $$

# This function involves a sum of terms, where each term depends on two adjacent variables. 
# The aim of optimization is to find the global minimum of this function within an $n$-dimensional space. 
# The Dixon-Price function is notable for having a global minimum at the point where all variables are equal to one.

# However, as the dimensionality increases beyond two, the Dixon-Price function exhibits multiple local minima, 
# making it a challenging optimization problem. This characteristic provides a valuable test for assessing the efficacy
# and robustness of optimization algorithms across various applications.




###############  DIXON-PRICE FUNCTION (OBJECTIVE) ###############

function dix(x::Vector) # x is a vector of real coordinates
    n = length(x) # Vector dimension
    sum = (x[1] - 1)^2 # First part of the function 
    for i in 2:n  # For each element 'x[i]' of the vector 'x'...
        sum += i * (2*x[i]^2 - x[i-1])^2 # Part of the function involving the summation
    end
    return sum # Sum of the first part with the second part
end


############### DIXON-PRICE FUNCTION (GRADIENT) ###############

function dpg(x::Vector) # x is a vector of real coordinates
    n = length(x)   # Vector dimension
    g = zeros(n)    # Creates a vector of n coordinates initialized to 0
    g[1] = 2 * (x[1] - 1) - 4 * (2 * x[2]^2 - x[1])  # Update the first coordinate of the vector with this result
    for i in 2:(n-1)
        g[i] = 2 * i * (2 * x[i]^2 - x[i-1]) * (4 * x[i]) - (i+1) * 2 * (2 * x[i+1]^2 - x[i]) # Update coordinates from 2 to n-1 with these results 
    end
    g[n] = n * 2 * (2 * x[n]^2-x[n-1]) * (4 * x[n]) # Update the last coordinate with this result
    return g
end


############### DIXON-PRICE FUNCTION (HESSIAN) ###############

function dph(x::Vector)   # x is a vector of real coordinates
    n = length(x)   # Vector dimension
    H = zeros(length(x), length(x))  # Creates a matrix of zeros with dimensions n x n 

    # For the first row of the matrix
    H[1, 1] = 6  # Replace the first result of the row with 6
    H[1, 2] = -16 * x[2] # Replace the second result of the row with this value

    # For rows 2 to n-1
    for i in 2:length(x)-1
        H[i, i-1] = -8 * i * x[i] # Replace the result of row i column i-1 with this value 
        H[i, i] = 8 * i * (2 * x[i]^2 - x[i-1]) + 32 * i * x[i]^2 + 2 * (i+1) # Replace the result of row i column i with this value
        H[i, i+1] = -8 * (i+1) * (x[i+1]) # Replace the result of row i column i+1 with this value
    end

    # For the last row
    H[end, end-1] = -8 * n * x[end]  # Replace the result of row n column n-1 with this value
    H[end, end] = 8 * n * (6 * x[end]^2 - x[end-1]) # Replace the result of row n column n with this value

    return H
end






############### RASTRING FUNCTION - DESCRIPTION ###############

# The Rastrigin function is a widely used multimodal minimization problem in optimization testing. 
# It is characterized by having numerous local minima that are uniformly distributed throughout its search space, 
# presenting a significant challenge for optimization algorithms.
# The Rastrigin function is defined as follows for $n$ dimensions:

# $$ f(x) = 10 \cdot n + \sum_{i=1}^{n} [x_i^2 - 10 \cdot \cos(2\pi x_i)] $$

# The function is typically evaluated within the hypercube $x_i \in [-5.12, 5.12]$ for all $i = 1, \ldots, n$. 
# These bounds define the range of each dimension within the search space.

# Despite its complexity, the global minimum of the Rastrigin function is known to occur at the origin $(0, 0, \ldots, 0)$, 
# where the function value is zero. However, finding this global minimum is challenging due to the presence of numerous
# local minima scattered throughout the search space.




############### RASTRING FUNCTION (OBJECTIVE) ###############

function rt(x::Vector{T}) where T <: Real  # x is a vector of real coordinates
    n = length(x)   # Vector dimension
    sum = 0  # First part of the function
    for i in 1:length(x)  # For each element 'x[i]' of the vector 'x'...
        sum += x[i]^2 -10*cos(2*pi*x[i])  # Part of the function involving the summation
    end
    a = 10 * n  # Part outside the summation
    return sum + a # Sum of the parts
end


############### RASTRING FUNCTION (GRADIENT) ###############
function rtg(x::Vector{T}) where T <: Real  # x is a vector of real coordinates
    g = zeros(length(x))  # Creates a vector of n coordinates initialized to 0
    n = length(x)  # Vector dimension
    for i in 1:length(x) # For each element 'x[i]' of the vector 'x'...
        g[i] = 2 * x[i] + 10 * sin(2π * x[i]) * 2π   # Update coordinates from 1 to n with these results 
    end
    return g 
end


############### RASTRING FUNCTION (HESSIAN) ###############

function rth(x::Vector)  # x is a vector of real coordinates
    n = length(x) # Vector dimension
    H = zeros(length(x), length(x)) # Creates a matrix of zeros with dimensions n x n 
    for i in 1:length(x)
        H[i, i] = 2 * (1 + 20 * π^2 * cos(2π * x[i]))  # Replace the result of row i column i with this value
    end
    return H
end






############### SPHERE FUNCTION - DESCRIPTION ###############

# The Sphere function is a fundamental optimization benchmark commonly used to evaluate optimization algorithms. 
# In $n$ dimensions, the Sphere function is defined as:

# $$ f(x) = \sum_{i=1}^{n} x_i^2 $$

# This function calculates the sum of squares of each variable, aiming to minimize the sum to achieve the global minimum. 
# The global minimum of the Sphere function is at the origin, where all variables are zero. The function represents a convex
# and smooth surface without local minima, making it relatively easier to optimize compared to many other functions.




############### SPHERE FUNCTION (OBJECTIVE) ###############

function sph(x::Vector{T}) where T <: Real # x is a vector of real coordinates
    n = length(x) # Vector dimension
    sum = 0   # First part of the function
    for i in 1:length(x)  # For each element 'x[i]' of the vector 'x'...
    sum += (x[i]^2)  # Part of the function involving the summation
    end
    return sum # Sum of the parts
end


############### SPHERE FUNCTION (GRADIENT) ###############

function sphg(x::Vector{T}) where T <: Real # x is a vector of real coordinates
    g = zeros(length(x)) # Creates a vector of n coordinates initialized to 0
    for i in 1:length(x) # For each element 'x[i]' of the vector 'x'...
        g[i] = 2 * x[i]  #Update coordinates from 1 to n with these results
    end
    return g 
end


############### SPHERE FUNCTION (HESSIAN) ###############

function sphh(x::Vector{T}) where T <: Real # x is a vector of real coordinates
    n = length(x)  # Vector dimension
    H = zeros(length(x), length(x))   # Creates a matrix of zeros with dimensions n x n 
    for i in 1:length(x)
        H[i, i] = 2    # Replace the result of row i column i with this value
    end
    return H
end






############### TRID FUNCTION - DESCRIPTION ###############

# The Trid function serves as a significant benchmark for evaluating optimization algorithms. Defined by the expression:

# $$ f(x) = \sum_{i=1}^{n} (x_i - 1)^2 - \sum_{i=2}^{n} x_i x_{i-1} $$

# This function calculates the sum of squares of each variable minus the sum of products of adjacent variables, 
# aiming to minimize this sum to achieve the global minimum. The global minimum of the Trid function occurs when 
# all variables are equal to their respective indices, resulting in a minimum value of 0. The function is convex and smooth, 
# devoid of local minima, which makes it relatively easier to optimize compared to many other functions.




############### TRID FUNCTION (OBJECTIVE) ###############

function tri(x::Vector{T}) where T <: Real  # x is a vector of real coordinates
    s1 = 0   # Initialize the variable 's1' to 0. 
    for i in 1:length(x)  # For each element 'x[i]' of the vector 'x'...
        s1 += (x[i] - 1)^2  # ...subtract 1 from 'x[i]', square the result, and add it to 's1'.
    end

    s2 = 0   # Initialize the variable 's2' to 0.
    for i in 2:length(x)  # For each element 'x[i]' of the vector 'x', starting from the second element...
        s2 += x[i] * x[i-1]  # ...multiply 'x[i]' by the previous element 'x[i-1]' and add the result to 's2'.
    end
   
    return s1 - s2  # Subtract 's2' from 's1' and return the result. This is the value of the 'trid' function for the input vector 'x'.
end





############### TRID FUNCTION (GRADIENT) ###############

function tridg(x::Vector{T}) where T <: Real
    n = length(x)   # Vector dimension
    g = zeros(n)    # Creates a vector of n coordinates initialized to 0
    g[1] = 2 * (x[1] - 1) - x[2]  # Update the first coordinate of the vector with this result
    for i in 2:length(x)-1 # For each element 'x[i]' of the vector 'x'...
        g[i] = 2 * (x[i] - 1) - (x[i-1] + x[i+1])  # Update coordinates from 2 to n-1 with these results 
    end
    g[end] = 2 * (x[n]-1) - x[n-1]
    return g
end


############### TRID FUNCTION (HESSIAN) ###############

function tridh(x::Vector)   # x is a vector of real coordinates
    n = length(x)   # Vector dimension
    H = zeros(length(x), length(x))  # Creates a matrix of zeros with dimensions n x n 

    # For the first row of the matrix
    H[1, 1] = 2  # Replace the first result of the row
    H[1, 2] = -1 # Replace the second result of the row with this value

    # For rows 2 to n-1
    for i in 2:length(x)-1
        H[i, i-1] = -1 # Replace the result of row i column i-1 with this value 
        H[i, i] = 2 # Replace the result of row i column i with this value
        H[i, i+1] = -1 # Replace the result of row i column i+1 with this value
    end

    # For the last row
    H[end, end-1] = -1  # Replace the result of row n column n-1 with this value
    H[end, end] = 2 # Replace the result of row n column n with this value

    return H
end 






############### SUM SQUARES FUNCTION - DESCRIPTION ###############

# The Sum Squares function stands as a fundamental benchmark used to assess optimization algorithms. 
# In $n$ dimensions, the function is expressed as:

# $$ f(x) = \sum_{i=1}^{d} i x_i^2 $$

# This function calculates a weighted sum of squares for each variable, with the weight increasing linearly with the index i. 
# The primary objective is to minimize this sum to attain the global minimum. The global minimum of the Sum Squares function is 
# achieved when all variables are zero, leading to a minimum value of 0. Notably, the function is characterized by its convex 
# and smooth nature, devoid of any local minima, which renders optimization relatively straightforward compared to other functions.




############### SUM SQUARES FUNCTION (OBJECTIVE) ###############

function squa(x::Vector{T}) where T <: Real
    n = length(x) # Vector dimension
    sum = 0   # First part of the function
    for i in 1:length(x)  # For each element 'x[i]' of the vector 'x'...
    sum += i * (x[i]^2)  # Part of the function involving the summation
    end
    return sum # Sum of the parts
end


############### SUM SQUARES FUNCTION (GRADIENT) ###############

function squagrad(x::Vector{T}) where T <: Real
    n = length(x)
    g = zeros(n) # Creates a vector of n coordinates initialized to 0
    for i in 1:n
        g[i] = 2 * i * x[i]  #Update coordinates from 1 to n with these results
    end
    return g 
end


############### SUM SQUARES FUNCTION (HESSIAN) ###############

function squahess(x::Vector{T}) where T <: Real
    n = length(x)  # Vector dimension
    H = zeros(n, n)   # Creates a matrix of zeros with dimensions n x n 
    for i in 1:n
        H[i, i] = 2 * i   # Replace the result of row i column i with this value
    end
    return H
end






############### STYBLINSKI-TANG FUNCTION - DESCRIPTION ###############

#
#
#
#



############### STYBLINSKI-TANG FUNCTION (OBJECTIVE) ###############

function stangobj(x::Vector{T}) where {T <: Real}
    # Inicializando a variável 'soma' com 0.0; ela irá armazenar a soma cumulativa dos cálculos
    soma = zero(T)
    
    # Iterando sobre cada índice do vetor x usando um loop for
    for i in 1:length(x)
        # Calculando o valor da função Styblinski-Tang para cada elemento x[i] em x e adicionando à 'soma'.
        # O cálculo envolve operações polinomiais em x[i].
        soma += (x[i]^4 - 16*x[i]^2 + 5*x[i])
    end
    
    # Multiplicando a soma final por 1/2 conforme a fórmula e retornando o resultado
    return soma / 2.0 
end




############### STYBLINSKI-TANG FUNCTION (GRADIENT) ###############

# Definindo a função gradiente da Styblinski-Tang
function stangrad(x::Vector{Float64})
    n = length(x)
    # Inicializando o vetor 'grad' com zeros; ele irá armazenar o gradiente
    grad = zeros(n)
    
    # Iterando sobre cada índice do vetor x usando um loop for
    for i in 1:n
        # Calculando o valor do gradiente da função Styblinski-Tang para cada elemento x[i] em x e armazenando no vetor 'grad'.
        # O cálculo envolve operações polinomiais em x[i].
        grad[i] = 1/2 * (4*x[i]^3 - 32*x[i] + 5)
    end
    
    # Retornando o vetor gradiente
    return grad
end



############### STYBLINSKI-TANG FUNCTION (HESSIAN) ###############

function stanghess(x::Vector{Float64})
    n = length(x)
    # Inicializando a matriz 'H' com zeros; ela irá armazenar a matriz Hessiana
    H = zeros(n, n)
    
    # Iterando sobre cada índice do vetor x usando um loop for
    for i in 1:n
        # A segunda derivada da função Styblinski-Tang em relação a x[i] é 1/2 * (12*x[i]^2 - 32).
        H[i, i] = 1/2 * (12*x[i]^2 - 32)
    end
    
    # Retornando a matriz Hessiana
    return H
end














# How to call the functions
stang = otf(stangobj, stangrad, stanghess)
sphere = otf(sph, sphg, sphh)
rast  = otf(rt, rtg, rth)
dp  = otf(dix, dpg, dph)
rosen = otf(rf, rg, rh)
trid = otf(tri, tridg, tridh)
squares = otf(squa, squagrad, squahess)

end