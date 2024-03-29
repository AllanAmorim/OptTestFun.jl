module OptTestFun

export quad1, quad2, rosen


# Write your package code here...
struct st_funlib
    obj  :: Function
    grad :: Function
    hess :: Function 
end


############ FUNCTION LIBRARY ###############################

function rf(x::Vector{T}) where T <: Real
    sum = 0 
    for i in 1:length(x)-1 
        sum += 100 * (x[i+1] - x[i]^2)^2 + (1 - x[i])^2 
    end
    return sum 
end

function rg(x::Vector{T}) where T <: Real
    n = length(x)
    g = zeros(n)
    g[1] = 200 * (x[2] - x[1]^2) * (-2 * x[1]) + 2 * (x[1] - 1)
    for i in 2:length(x)-1
        g[i] = 200 * (x[i] - x[i-1]^2) - 400 * x[i] * (x[i+1] - x[i]^2) + 2 * (x[i] - 1)
    end
    g[end] = 200.0 * (x[end] - x[end-1]^2)
    return g
end

function rh(x::Vector{T}) where T <: Real
    n = length(x)
    H = zeros(length(x), length(x)) 
    H[1, 1] = -400 * (x[2] - x[1]^2) + 800 * x[1]^2 + 2
    H[1, 2] = -400.0 * x[1] 
    for i in 2:length(x)-1 
        H[i, i-1] = -400 * x[i-1]
        H[i, i] = 200 - 400 * (x[i+1] - x[i]^2) + 800 * x[i]^2 + 2
        H[i, i+1] = -400.0 * x[i] 
    end
    H[end, end-1] = -400.0 * x[end-1] 
    H[end, end] = 200.0 
    return H  
end


function f1(x)
    return sum(x)
end

function g1(x)
    return 2x
end

function h1(x)
    n = size(x,1)
    h = zeros(n,n)
    for i in 1:n
        h[i,i] = 2
    end
    return h
end


function f2(x)
    return 3*sum(x)
end

function g2(x)
    return 5x
end

function h2(x)
    n = size(x,1)
    h = zeros(n,n)
    for i in 1:n
        h[i,i] = 2*x[i]
    end
    return h
end

# Function 

rosen = st_funlib(rf, rg, rh)
quad1 = st_funlib(f1,g1,h1)
quad2 = st_funlib(f2,g2,h2)

end