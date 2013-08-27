function elev_goose, ioff, phidiff, maxbeam, bmsep, bmnum, tfreq, tdiff, phi0,$
										debug=debug

C = 2.997e8   ; speed of light in free space

; ioff[3] : interferometer offset [m]
;           +x along array towards increasing beam numbers
;           +y boresight direction with sign determined by x and z
;           +z up

;  The coordinate system is: +x is in front of the main array,
;  +y is to the left of the center, and +z is upward.   

;  The angle xi is in the x-y plane.  The angle is measured counter-clockwise
;  from the x axis.  If the y offset and x offset are both positive, the
;  angle is positive.  

;  Sep is the distance between the centers of the two arrays

;  If the z offset is not zero, then the elevation angle has to be corrected.
;  An apparent elevation angle of 0 would actually being pointing slightly
;  upward (downward) if the z offset is positive (negative).

xi        = atan2(ioff[0],ioff[1])
sep       = sqrt(ioff[0]*ioff[0] + ioff[1]*ioff[1] + ioff[2]*ioff[2])

; note that this correction is atan(z/sep) in elevation.pro
elev_corr = atan2(ioff[2],sqrt(ioff[0]*ioff[0] + ioff[1]*ioff[1]))

cos_xi    = cos(xi)
sin_xi    = sin(xi)
cos_xi2   = cos_xi*cos_xi
sin_xi2   = sin_xi*sin_xi
  
; compute phasing matrix cone angle
offset = maxbeam/2. - 0.5
psi    = bmsep*(bmnum - offset)*!pi/180.	; note, this is phi in elevation.pro

; compute wavenumber
k = 2. * !pi * tfreq * 1000.0/C

; delay due to path difference
dchi_cable = - 2 * !pi * tfreq * 1000. * tdiff * 1.e-6

; compute the minimum cone angle (alpha)

temp          = sin(psi) + sin_xi
sin_psi_xi    = sin(psi)*sin_xi

sn2_alpha_min = (temp*temp)/(1. + sin_psi_xi)
cs_alpha_min  = sqrt(1. - sn2_alpha_min);

; now iterate sn2_alpha_min, cs_alpha_min to improve value of alpha_min

sin_psi2      = sin(psi)*sin(psi)
sin_psi2_xi2  = sin_psi2 + sin_xi2

sn2_alpha_old = sn2_alpha_min

sn2_alpha_min = sin_psi2_xi2 + 2.0*cs_alpha_min*sin_psi_xi
while (abs(sn2_alpha_min - sn2_alpha_old) gt 0.005*sn2_alpha_old) do begin
	cs_alpha_min = sqrt(1.0 - sn2_alpha_min)
	sn2_alpha_old = sn2_alpha_min
	sn2_alpha_min = sin_psi2_xi2 + 2.0*cs_alpha_min*sin_psi_xi
endwhile
;repeat begin
;	sn2_alpha_min = sin_psi2_xi2 + 2.0*cs_alpha_min*sin_psi_xi
;	cs_alpha_min  = sqrt(1. - sn2_alpha_min)
;	sn2_alpha_old = sn2_alpha_min
;	if keyword_set(debug) then print, cs_alpha_min, sn2_alpha_min
;stop
;end until (abs(sn2_alpha_min - sn2_alpha_old) gt 0.005*sn2_alpha_old)

cs_alpha_min  = sqrt(1. - sn2_alpha_min)

; we've now got the sin & cos of alpha_min
; compute the total phase difference

dchi_sep_max  = k * sep / cos_xi * cs_alpha_min
dchi_max      = dchi_cable + dchi_sep_max
;n             = 0.5 - dchi_max/(2*PI)
; this needs to be an integer, but it isn't clear what would happen in C
;n             = round(0.5 - dchi_max/(2*!pi))
n             = fix(0.5 - dchi_max/(2*!pi))

; this should be the true phase difference
dchi          = phi0 - n*2*!pi

if (dchi gt dchi_max) then             dchi -= 2*!pi
if (dchi lt (dchi_max - (2*!pi))) then dchi += 2*!pi

; compute the cone angle (alpha)

dchi_old = 0.
sn2_eps  = 0.
while (abs(dchi_old - dchi) gt !pi) do begin
	cs_alpha = (dchi - dchi_cable)/(k*sep)*cos_xi
	sn2_eps  = 1. - (cs_alpha*cs_alpha)/(cos_xi2) - (sin_psi2/cos_xi2) $
								- 2.*cs_alpha*sin_psi_xi/cos_xi2
	dchi_old = dchi

	if ((abs(sn2_eps) gt 1.) or (sn2_eps lt 0.)) then begin
		dchi -= 2*!pi
	  print, 'changing dchi by -2pi. ',dchi_old,' -> ',dchi
		stop
	endif
endwhile

sn_eps = sqrt(sn2_eps)
elev   = asin(sn_eps)

; The calculated elevation angle is actually with respect to the plane
;   that includes the two antenna arrays. This has to be corrected for the
;   difference in elevation between the front array and the back array.

elev = elev + elev_corr
return, 180.*elev/!pi

end

