; test for Ray
m2ft    = 3.28084	; feet per meter
theta   = 0.			; elevation angle
bmnum   = 11
tfreq   = 10000.	; 10 MHz
ioff    = [0., -80., 0]		; CV offsets
tdiff   = -0.351	; in us

C       = 2.997e8
maxbeam = 24
bmsep   = 3.24
offset  = maxbeam/2. - 0.5
phi     = bmsep * (bmnum - offset) * !pi/180.

k = 2*!pi * tfreq * 1000./C		; wave number in free space

cphi = cos(phi)
sphi = sin(phi)
cthe = cos(theta*!pi/180.)
sthe = sin(theta*!pi/180.)

A = cthe*sphi
B = cthe*cphi
C = sthe

d_i = A*ioff[0] + B*ioff[1] + C*ioff[2]
d_m = 0 ; always make it at the origin....

d = d_i-d_m

print, 'Compute Phi0'
print, '------------'
print, 'Phi        = ', phi*180./!pi
print, 'separation = ', d, ' m'
print, '           = ', d*m2ft, ' ft'

phi_geo = k*d
print, 'Phi_geo    = ', phi_geo, ' rad'

phi_path = 2*!pi*tfreq*1e3*tdiff*1e-6
print, 'Phi_path   = ', phi_path, ' rad'

Phi_total = phi_path - phi_geo
print, 'Phi_total  = ', phi_total, ' rad'


end

