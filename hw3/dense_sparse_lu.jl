using LinearAlgebra
using Plots
using Printf # for formatting text output

"""
    var_set(N)

Initializes banded sparse matrix, and vector b of size N


A = -2 1 0 .... 0
     0-2 1 0 .. 0
     0 0-2 1 0  0
"""
function var_set(N)
    A = Matrix{Float64}(undef, N, N)
    A = zeros(N,N)

    A_d = Matrix{Float64}(undef, N, N)
    A_d = ones(N,N)
    for i = 1:N
        A_d[i,i] = -2
    end

    b = rand(N,1)
    b = b[:]

    for i = 1:N
        b[i] = 1
    end

    for i = 1:N
        A[i,i] = -2
        if i != N
            A[i, i+1] = 1
            A[i+1, i] = 1
        end
    end

    return (A, A_d, b)
end


"""
    LU_solve(A, b)

Uses Julia linear algebra package to solve Ax = b
given Matrix A and vector b
"""
function LU_solve(A, b)
    F = lu(A)
    # display(F)

    x = F\b

    # display(x)
    return x
end

function main()
    dense_time = []
    sparse_time = []
    testSizes = [10, 100, 1000, 10000]
    @printf("starting loop for values")
    for i = 1:size(testSizes,1)
        (A, A_d, b) = var_set(testSizes[i])
        temp = @timed x = LU_solve(A, b)
        push!(sparse_time, temp[2])
        new_temp = @timed x = LU_solve(A_d, b)
        push!(dense_time, new_temp[2])
        @printf("done at N = %d\n", testSizes[i])
    end
    scatter(testSizes, dense_time, color="red", xlabel="Size of N", ylabel="Time", title = "Size of N vs Time", legend=false)
    scatter!(testSizes, sparse_time, color="blue", legend=false)
    savefig("denseSparsePlot.png")
end

main()
