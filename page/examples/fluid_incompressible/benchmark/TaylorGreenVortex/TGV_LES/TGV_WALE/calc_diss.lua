-- declear data array
nVals = 0
ke   = {}
diss = {}
time = {}

-- get KE results file
filename = os.getenv("ke_file")
print("Open file: "..filename)
local f = assert(io.open(filename, "r" ))

-- read in the first two lines which are useless
print("Read the first two lines.")
t = f:read("*line")
t = f:read("*line")

-- start read all data
while true do
  -- first 3 columns are coordinates
  local timeVal, keVal = f:read("*number", "*number")

  -- check if it already reaches the end of file (baryX == nil)
  if not timeVal then break end

  -- count number of values
  nVals = nVals + 1
  -- append val to array
  ke[nVals] = keVal
  time[nVals] = timeVal
end
-- close file
f:close()

-- calculate -dE/dt = ( E(t) - E(t+1) ) / dt
print("The number of data = "..nVals)

filename = os.getenv("dr_file")
f = assert(io.open(filename, "w" ))
for i = 1,nVals-1 do
  diss[i] = (ke[i] - ke[i+1]) / (time[i+1] - time[i])
  f:write(string.format("%16.11f",   time[i] ))
  f:write(string.format("%16.11f\n", diss[i] ))
end
-- close file
f:close()
