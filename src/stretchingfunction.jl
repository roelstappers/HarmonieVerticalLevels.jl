using Zygote, Plots, DelimitedFiles, Dierckx, TOML

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


# LLAPRXPK:
# Full layers are assumed to be computed as for the options LVERTFE=F, NDLNPR=0 of ARPEGE/ALADIN.
# LLAPRXPK=T => pressure(l)=0.5(pressure(lbar-1)+pressure(lbar))
#   ("l" stands for full levels, "lbar" for half levels).
# LLAPRXPK=F => a more tricky way to compute pressure(l).
# When using the vertical layers for LVERTFE=F, NDLNPR=0, LAPRXPK=F
# in the model, it is recommended to use LLAPRXPK=F.
# When using the vertical layers for LVERTFE=F, NDLNPR=0, LAPRXPK=T
# of for LVERTFE=T, it is recommended to use LLAPRXPK=T.
const LLAPRXPK = true 
δπ₁ = LLAPRXPK ? 2π₁ : ℯ * π₁ 


#  α₁,α₃,αₕ: coefficients defining the function allowing
#  to compute the A and B.
#  For α₁ and α₃ it is recommended to take values between 1 and 5.
#  For αₕ it is recommended to take values between -3 and -1
# (αₕ must bever be > -1).

α₁ =  2.8; α₃ =  1.7

x₁ = 1.0 / N; 
x₂ = NSTR / N
x₃ = (N - NPBL) / N
x₄ = (N - 1.0) / N

y₁ = δπ₁ / π₀₀
y₂ = πₛₜᵣ / π₀₀ 
y₃ = πPBL / π₀₀
y₄ = (π₀₀ - δπL) / π₀₀


d₁ = (x₁ * y₂ - x₂ * y₁) * (1 / x₁) * (x₂ - x₁)^(-α₁) 
d₃ = ((1 - x₄) * (1 - y₃) - (1 - x₃) * (1 - y₄)) / (1 - x₄) * (x₄ - x₃)^(-α₃)


"""
Definition of the Stretching function. 
See Pierre Bénard section 3 
"""
function m(x)  
    ystr(x) = (y₁ / x₁) * x + d₁ * (x - x₁)^α₁         
    ypbl(x) = 1 - (1 - y₄) / (1 - x₄) * (1 - x) - d₃ * (x₄ - x)^(α₃)    
    if 0.0 ≤ x ≤ x₁      # Upper most sub-domain
        return y₁ / x₁ * x
    elseif x₁ < x ≤  x₂  # "Strato" sub-domain               
        return ystr(x)
    elseif x₂ < x ≤ x₃   # "Tropo sub-domain        
        Δx = x₃ - x₂; Δy = y₃ - y₂; s = Δy / Δx
        return y₂ + (x - x₂) * ystr'(x₂)  + (x - x₂)^2 * (Δx * (s - ystr'(x₂)) + (x - x₃) * (ystr'(x₂) + ypbl'(x₃) - 2s)) / (Δx^2)
    elseif x₃ < x ≤ x₄   # "PBL" sub-domain        
        return ypbl(x) 
    elseif x₄ < x ≤ 1.0   # Bottom sub domain
        return 1.0 - (1.0 - y₄) / (1.0 - x₄) * (1 - x)
    else 
        DomainError(x, "m is only defined for 0 ≤ x ≤ 1")
    end 
end



αₕ = -1.6

"""
Definition of the hybridicity function. 
See Pierre Bénard section 4 
"""
function h(y)
    yσ = m((N - Nσ) / N)
    yπ = m(Nπ / N) 
    d₁ = αₕ * yσ^2 / (yσ - yπ)   
    d₂ = 1 + αₕ * yσ / (yσ - yπ)  
    if 0.0 ≤ y ≤ yπ 
       0.0 
    elseif yπ ≤ y ≤ yσ
       d₁ / (d₂ - ((y - yπ) / (yσ - yπ))^αₕ) 
    elseif yσ ≤ y ≤ 1.0
       y      
    else
        DomainError(y, "h is only defined for 0 ≤ y ≤ 1")
    end    
end 

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