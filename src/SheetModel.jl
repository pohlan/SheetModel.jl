module SheetModel

alpha = 1.25
beta = 1.5
function calc_q(h, dphi_du, k) # u can be x or y
    q = - k .* h.^alpha .* abs.(dphi_du).^(beta-2) .* dphi_du
    return q
end

end # module
