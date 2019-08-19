using Zygote, Plots, DelimitedFiles, Dierckx, TOML, LaTeXStrings

const N    = 65      # Total number of layers.
const NPBL = 8       # Total number of layers in the PBL.
const NSTR = 27      # Total number of layers in the stratosphere.
const Nπ   = 12      # Total number of pure pressure layers (minimum 1 included).
const Nσ   = 2       # Total number of pure sigma layers (minimum 1 included).

const πPBL = 90000.  # pressure of the top of the PBL in pascals
const πₛₜᵣ = 25000.  # pressure of the tropopause in pascals
const π₁   = 9.9     # pressure of the full layer l=1, in pascals
const δπL  = 410.    # [Delta pressure] of the full layer l=L, in pascals 
const π₀₀  = 101325. # standard surface pressure, in pascals




#  α₁,α₃,αₕ: coefficients defining the function allowing
#  to compute the A and B.
#  For α₁ and α₃ it is recommended to take values between 1 and 5.
#  For αₕ it is recommended to take values between -3 and -1
# (αₕ must bever be > -1).

α₁ =  2.8; α₃ =  1.7





αₕ = -1.6

A(x) = π₀₀ * (m(x) - h(m(x)))   # Bénard equation  14
B(x) = h(m(x))                  # Bénard equation  15



x = range(0.0, 1.0, length = N+1)
writedlm("A.txt", A.(x))
writedlm("B.txt", B.(x))

# plot(A',0,1)

plot(m, 0, 1, label = "Bénard")
scatter!([x₁, x₂, x₃, x₄], [y₁, y₂, y₃, y₄], label = "")
plot!(yticks = ([y₁, y₂, y₃, y₄],  [L"y_1", L"y_2", L"y_3",L"y_4"]))
plot!(xticks = ([x₁, x₂, x₃, x₄],  [L"x_1", L"x_2", L"x_3",L"x_4"]))
plot!(legend = :topleft)
plot!(title = "Stretching function")

# Compare with cubic splines 

x = [0.0, x₁, x₂, x₃, x₄, 1.0]
y = [0.0, y₁, y₂, y₃, y₄, 1.0]
sp1 = Spline1D(x, y)
plot!(x->sp1(x), 0, 1, label = "Spline")


#plot hybridicity function
yσ = m((N - Nσ) / N)
yπ = m(Nπ / N) 
plot(h, 0, 1) 
plot!(x->x, 0, 1, linestyle = :dash)
plot!([yπ  yσ], [0 yσ], seriestype = :scatter)
plot!(xticks = ([0, yπ,  yσ, 1], [L"0",L"y\pi",L"y\sigma",L"1"]))
plot!(yticks = ([0, 1],[L"0",L"1"]) )
plot!(legend = false)


#plot(A,0,1)
#VLEV = TOML.parsefile("../Harmonie.jl/test/config/VLEV/65.toml")
#VLEV = TOML.parsefile("../Harmonie.jl/test/config/VLEV/ECMWF_60.toml")
#plot!(range(0,1,length=VLEV["NLEV"]+1), VLEV["AHALF"])

# end