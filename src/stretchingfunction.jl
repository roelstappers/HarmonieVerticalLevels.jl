
# module Verticallevels
using Zygote, Plots 

const N     = 55 # Total number of layers.
const NPBL  = 8  # Total number of layers in the PBL.
const NSTR  = 27 # Total number of layers in the stratosphere.
const NPRES = 12 # Total number of pure pressure layers (minimum 1 included).
const NSIGM = 1  # Total number of pure sigma layers (minimum 1 included).

const πPBL = 90000. # pressure of the top of the PBL in pascals
const πₛₜᵣ = 25000. # pressure of the tropopause in pascals
const π₁  = 9.9         # pressure of the full layer l=1, in pascals
const δπL  = 410. # [Delta pressure] of the full layer l=L, in pascals 
const π₀₀  = 101325.    # standard surface pressure, in pascals


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

x₁=1.0/N; 
x₂=NSTR/N
x₃=(N-NPBL)/N
x₄=(N-1.0)/N

y₁=δπ₁/π₀₀
y₂=πₛₜᵣ/π₀₀ 
y₃=πPBL/π₀₀
y₄=(π₀₀-δπL)/π₀₀


d₁ = (x₁*y₂-x₂*y₁)*(1/x₁)*(x₂-x₁)^(-α₁) 
d₃ = ((1-x₄)*(1-y₃) - (1-x₃)*(1-y₄))/(1-x₄)*(x₄-x₃)^(-α₃)

ystr(x) = (y₁/x₁)*x + d₁*(x-x₁)^α₁         
ypbl(x) = 1-(1-y₄)/(1-x₄)*(1-x)-d₃*(x₄-x)^(α₃)


function m(x)  
    if 0.0 ≤ x ≤ x₁      # Upper most sub-domain
        return y₁/x₁*x
    elseif x₁ < x ≤  x₂  # "Strato" sub-domain               
        return ystr(x)
    elseif x₂ < x ≤ x₃   # "Tropo sub-domain        
        Δx = x₃ - x₂; Δy = y₃ - y₂; s = Δy/Δx
        return y₂ + (x-x₂)*ystr'(x₂)  + (x- x₂)^2*(Δx*(s-ystr'(x₂)) + (x-x₃)*(ystr'(x₂) + ypbl'(x₃) -2s))/(Δx^2)
    elseif x₃ < x ≤ x₄   # "PBL" sub-domain        
        return ypbl(x) 
    elseif x₄ < x ≤ 1.0   # Bottom sub domain
        return 1.0 - (1.0-y₄)/(1.0-x₄)*(1-x)
    else 
        DomainError(x,"m is only defined for 0 ≤ x ≤ 1")
    end 
end

plot(m,0,1,label="Bénard")
scatter!([x₁, x₂, x₃, x₄],[y₁, y₂, y₃, y₄],label="")
plot!(yticks = ([y₁, y₂, y₃, y₄],  ["y1", "y2", "y3","y4"]))
plot!(xticks = ([x₁, x₂, x₃, x₄],  ["x1", "x2", "x3","x4"]))
plot!(legend=:topleft)
plot!(title="Stretching function")

# Compare with cubic splines 
using Dierckx
x = [0.0, x₁, x₂, x₃, x₄, 1.0]
y = [0.0, y₁, y₂, y₃, y₄, 1.0]
sp1 = Spline1D(x,y)
plot!(x->sp1(x),0,1,label="Spline")

# end