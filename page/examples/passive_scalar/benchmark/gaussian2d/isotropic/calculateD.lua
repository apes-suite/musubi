-- Calculate diffusion factor from Gaussian Hill profile at a certain
-- time step.

require 'args'
if order == 'first' then
    order_to_write = '1st'
else
    order_to_write = '2nd'
end
case_args = tostring(tau)..'_'..tostring(u_field)..'_'..order_to_write

case_name = 'spc1'
sim_name = 'simulation_'..case_name..'_p00000_t200.000E+00.res'
local file = io.open('tracking/'..sim_name)

local data = {}
local offset, amplititude = c_add, c0

for line in file:lines() do
    local column1, column2, column3, column4, column5 =
        line:match("(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s+(%S+)")

    column1 = tonumber(column1)
    column4 = tonumber(column4)
    column5 = tonumber(column5)

    if type(column1) == 'number' then
        local row = {
            column1,
            (column4 - offset) / amplititude / (2*math.pi) / sigma0^2,
            (column5 - offset) / amplititude / (2*math.pi) / sigma0^2
        }
        table.insert(data, row)
    end
end

file:close()

-- for _, row in ipairs(data) do
--     print(table.concat(row, "  "))
-- end

local aves = {D_sim = 0., D_real = 0.}

for _, row in ipairs(data) do
    aves.D_sim = aves.D_sim + row[2] * row[1]
    aves.D_real = aves.D_real + row[3] * row[1]
end

local sigmas = {D_sim = 0., D_real = 0.}

for _, row in ipairs(data) do
    sigmas.D_sim = sigmas.D_sim + (row[1] - aves.D_sim)^2 * row[2]
    sigmas.D_real = sigmas.D_real + (row[1] - aves.D_real)^2 * row[3]
end

D = {
    sim = (sigmas.D_sim - sigma0 ^ 2) / 2 / t_total,
    real = (sigmas.D_real - sigma0 ^ 2) / 2 / t_total,
}

print(
    string.format('%-20s', case_args)..'\t'..
    string.format('%20.18f', tau)..'\t'..
    string.format('%20.18f', u_field)..'\t'..
    string.format('%20.18f', D.sim)..'\t'..
    string.format('%20.18f', D.real)..'\t'..
    string.format('%20.18f', (D.sim - D.real) / D.real)..'\t'..
    string.format('%20.18f', (D.sim - D.real))
)

