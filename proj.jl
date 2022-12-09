using Plots
using CUDA

const n = 2^20
const THREADS_PER_BLOCK = 256

C = zeros(Float32, n) |> cu
A = fill(1.0f0, n) |> cu
B = fill(2.0f0, n) |> cu

function add!(c, a, b)
    x = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    @inbounds c[x] = a[x] + b[x]
    return
end

threads = THREADS_PER_BLOCK # 256
blocks = n ÷ threads # 4096

# launch the kernel with 4096 blocks, of 256 threads each
print("starting")
@cuda threads=threads blocks=blocks add!(C, A, B)
print("done")
