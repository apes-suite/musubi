require 'params'

-- Coefficients for the approximation
local a1 =  0.254829592
local a2 = -0.284496736
local a3 =  1.421413741
local a4 = -1.453152027
local a5 =  1.061405429
local p  =  0.3275911

-- Sign function
local function sign(x)
    return x < 0 and -1 or 1
end

-- Approximation of erfc(x)
local function erfc(x)
    -- Save the sign of x
    local sign_x = sign(x)
    x = math.abs(x)
    
    -- Calculate t
    local t = 1 / (1 + p * x)
    
    -- Use Horner's method to calculate the approximation
    local y = (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t
    
    -- Compute erfc
    local erfc_approx = y * math.exp(-x * x)
    
    -- Adjust for the sign of the input
    if sign_x == -1 then
        return 2 - erfc_approx
    else
        return erfc_approx
    end
end

-- Analytical solution by Bear(1978)
function bearSol(x, y, z, t)
    local t_mol = t + t_ini

    -- Compute beta
    local beta = math.sqrt(v^2 / 4 / D^2 + lambda / D)
    
    -- Compute x * beta
    local xbeta = x * beta

    -- Compute offset of efrc function
    local offset = math.sqrt(v^2 + 4 * lambda * D) * t_mol

    -- Compute denomitor of offset
    local deno = 2 * math.sqrt(D * t_mol)
    
    return C_0 / 2 * math.exp(v * x / 2 / D) * (math.exp(-xbeta) * erfc((x - offset) / deno) + math.exp(xbeta) * erfc((x + offset) / deno))
end
