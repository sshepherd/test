;-----------------------------------------------------------------------------
;
; function to compute elevation angle using separation between planes
; set 'corr' for beams direction as function of elevation angle
;

function elev_new, ioff,phidiff,maxbeam,bmsep,bmnum,tfreq,tdiff,phi_obs,$
									corr=corr

; compute beam direction, defined at zero elevation angle
offset = maxbeam/2. - 0.5
phi0   = bmsep * (bmnum - offset) * !dtor

X = ioff[0]
Y = ioff[1]
Z = ioff[2]

C = 2.997e8
k = 2*!pi * tfreq * 1000./C

; this term needs to change for more compicated geometries, i.e.,
;  maximum plane separation is no longer at zero elevation or phi0

; determine number of 2pis to shift phi_obs by
dmax = -(X*sin(phi0) + Y*cos(phi0))	; maximum distance between planes
																	; zero elevation angle

; phase delay [radians] due to electrical path difference
phi_path  = 2.*!pi * tfreq*1000. * tdiff*1e-6

phi_max = phi_path - dmax*k		; last term is max phi_geo

n2pi = fix((phi_obs-phi_max)/(2*!pi))		; flip for sign changes?
d2pi = n2pi*2*!pi
phi_obs -= d2pi

; NOTE: +/- phi_obs does not seem to matter, rather sign of dmax does...

; now solve for theta
if keyword_set(corr) then begin
	D = 1. - sin(phi0)*sin(phi0)
	E = (phi_path + phi_obs)/k - sin(phi0)*X
	sqrtfac = -E*E + (Z*Z+Y*Y)*D
	denom = Z*Z + Y*Y
	nump  = E*Z + abs(Y)*sqrt(sqrtfac)
	numn  = E*Z - abs(Y)*sqrt(sqrtfac)
endif else begin
	D = X*sin(phi0) + Y*cos(phi0)
	E = (phi_path + phi_obs)/k
	sqrtfac = D*D*(Z*Z + D*D - E*E)
	denom = Z*Z + D*D
	nump  = E*Z + sqrt(sqrtfac)
	numn  = E*Z - sqrt(sqrtfac)
endelse

theta = asin(nump/denom)/!dtor

return, theta
end

