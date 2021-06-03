using Pkg
Pkg.activate(@__DIR__)
using SheetModel
S = SheetModel

xc, yc, ϕ0, ϕ = S.runthemodel();
S.plot_output(xc, yc, ϕ0, ϕ)
