
using HarmonieVerticalLevels
using Plots

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