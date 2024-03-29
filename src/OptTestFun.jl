module OptTestFun

export quad1, quad2, rosen, dp, rast, sphere

struct st_funlib
    obj  :: Function
    grad :: Function
    hess :: Function 
end

############### ROSENBROCK FUNCTION ###############################

## OBJECTIVE FUNCTION
function rf(x::Vector{T}) where T <: Real
    # x is a vector of real coordinates
    n = length(x)   # Vector dimension
    sum = 0         # First part of the function
    for i in 1:length(x)-1 
        sum += 100 * (x[i+1] - x[i]^2)^2 + (1 - x[i])^2  # Part of the function involving the summation
    end
    return sum     # Sum of the first part with the second part
end

## GRADIENT FUNCTION
function rg(x::Vector{T}) where T <: Real
    # x is a vector of real coordinates
    n = length(x)   # Vector dimension
    g = zeros(n)    # Creates a vector of n coordinates initialized to 0
    g[1] = 200 * (x[2] - x[1]^2) * (-2 * x[1]) + 2 * (x[1] - 1)   # Update the first coordinate of the vector with this result
    for i in 2:length(x)-1
        g[i] = 200 * (x[i] - x[i-1]^2) - 400 * x[i] * (x[i+1] - x[i]^2) + 2 * (x[i] - 1) # Update coordinates from 2 to n-1 with these results 
    end
    g[end] = 200.0 * (x[end] - x[end-1]^2)  # Update the last coordinate with this result
    return g
end

## HESSIAN FUNCTION
function rh(x::Vector{T}) where T <: Real
    # x is a vector of real coordinates
    n = length(x)   # Vector dimension
    H = zeros(length(x), length(x))   # Creates a matrix of zeros with dimensions n x n

    # For the first row of the matrix
    H[1, 1] = -400 * (x[2] - x[1]^2) + 800 * x[1]^2 + 2  # Replace the first result in the line with this value
    H[1, 2] = -400.0 * x[1] # Replace the second result of the row with this value

    # For rows 2 to n-1
    for i in 2:length(x)-1 
        H[i, i-1] = -400 * x[i-1]   # Replace the result of row i column i-1 with this value
        H[i, i] = 200 - 400 * (x[i+1] - x[i]^2) + 800 * x[i]^2 + 2  # Replace the result of row i column i with this value
        H[i, i+1] = -400.0 * x[i]  # Replace the result of row i column i+1 with this value
    end

    # For the last row
    H[end, end-1] = -400.0 * x[end-1]  # Replace the result of row n column n-1 with this value
    H[end, end] = 200.0  # Replace the result of row n column n with this value

    return H  
end

############### DIXON-PRICE FUNCTION ###############################

## OBJECTIVE FUNCTION
function dix(x::Vector) # x is a vector of real coordinates
    n = length(x) # Vector dimension
    sum = (x[1] - 1)^2 # First part of the function 
    for i in 2:n
        sum += i * (2*x[i]^2 - x[i-1])^2 # Part of the function involving the summation
    end
    return sum # Sum of the first part with the second part
end

## GRADIENT FUNCTION
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

## HESSIAN FUNCTION
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

############### RASTRING FUNCTION ###############################

## OBJECTIVE FUNCTION
function rt(x::Vector{T}) where T <: Real  # x is a vector of real coordinates
    n = length(x)   # Vector dimension
    sum = 0  # First part of the function
    for i in 1:length(x)
        sum += x[i]^2 -10*cos(2*pi*x[i])  # Part of the function involving the summation
    end
    a = 10 * n  # Part outside the summation
    return sum + a # Sum of the parts
end

## GRADIENT FUNCTION
function rtg(x::Vector{T}) where T <: Real  # x is a vector of real coordinates
    g = zeros(length(x))  # Creates a vector of n coordinates initialized to 0
    n = length(x)  # Vector dimension
    for i in 1:length(x) 
        g[i] = 2 * x[i] + 10 * sin(2π * x[i]) * 2π   # Update coordinates from 1 to n with these results 
    end
    return g 
end

## HESSIAN FUNCTION
function rth(x::Vector)  # x is a vector of real coordinates
    n = length(x) # Vector dimension
    H = zeros(length(x), length(x)) # Creates a matrix of zeros with dimensions n x n 
    for i in 1:length(x)
        H[i, i] = 2 * (1 + 20 * π^2 * cos(2π * x[i]))  # Replace the result of row i column i with this value
    end
    return H
end

############### SPHERE FUNCTION ###############################

## OBJECTIVE FUNCTION
function sph(x::Vector{T}) where T <: Real # x is a vector of real coordinates
    n = length(x) # Vector dimension
    sum = 0   # First part of the function
    for i in 1:length(x) 
    sum += (x[i]^2)  # Part of the function involving the summation
    end
    return sum # Sum of the parts
end

## GRADIENT FUNCTION
function sphg(x::Vector{T}) where T <: Real # x is a vector of real coordinates
    g = zeros(length(x)) # Creates a vector of n coordinates initialized to 0
    for i in 1:length(x)
        g[i] = 2 * x[i]  #Update coordinates from 1 to n with these results
    end
    return g 
end

## HESSIAN FUNCTION
function sphh(x::Vector{T}) where T <: Real # x is a vector of real coordinates
    n = length(x)  # Vector dimension
    H = zeros(length(x), length(x))   # Creates a matrix of zeros with dimensions n x n 
    for i in 1:length(x)
        H[i, i] = 2    # Replace the result of row i column i with this value
    end
    return H
end

# How to call the functions
sphere = st_funlib(sph, sphg, sphh)
rast  = st_funlib(rt, rtg, rth)
dp  = st_funlib(dix, dpg, dph)
rosen = st_funlib(rf, rg, rh)

end
