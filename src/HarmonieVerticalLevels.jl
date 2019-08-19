

module HarmonieVerticalLevels

using Zygote, Plots, DelimitedFiles, Dierckx, TOML, LaTeXStrings


export h, m

include("stretchingfunction.jl")
include("hydricityfunction.jl")

A(x) = π₀₀ * (m(x) - h(m(x)))   # Bénard equation  14
B(x) = h(m(x))                  # Bénard equation  15



end