"""
Definition of the hybridicity function. 
See Pierre Bénard section 4 
"""
function h(y;yσ, yπ, αₕ)    
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
