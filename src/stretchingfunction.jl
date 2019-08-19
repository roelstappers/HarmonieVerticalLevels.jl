"""
Definition of the Stretching function. 
See Pierre Bénard section 3 
"""
function m(x; N,NSTR,NPBL, π₀₀, π₁, δπL,πPBL,πₛₜᵣ)

    const LLAPRXPK = true 
    δπ₁ = LLAPRXPK ? 2π₁ : ℯ * π₁ 

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
