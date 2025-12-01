theta = math.pi / 4.
n = 20


tau = 0.8
u_mag = 0.1
nelem = 100
length_usr = 256
level_usr = math.log(length_usr)/math.log(2)
t_total = 200
interval = t_total
sigma0 = 10
c_add = 0.1
c0 = (1.2-c_add)

u_field = {u_mag * math.cos(theta), u_mag * math.sin(theta)}
k_L = (tau - 0.5) / 3 / u_mag
k_T = k_L / 2^n

Dxx = u_mag * k_T + (k_L - k_T) * u_field[1]^2 / u_mag
Dyy = u_mag * k_T + (k_L - k_T) * u_field[2]^2 / u_mag
Dxy = (k_L - k_T) * u_field[1] * u_field[2] / u_mag
Dzz = u_mag * k_T

