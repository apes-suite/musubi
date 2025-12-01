require 'arg_given'

nelem = 100     -- cube length: 2*nelem
length_usr = 256
level_usr = math.log(length_usr)/math.log(2)
t_total = 200   -- maximum time
sigma0 = 10     -- initial sigma of gaussian hill
c_add = 0.1     -- minimum value of gaussian hill
c0 = (1.2-c_add)    -- amplitude of gaussian hill
