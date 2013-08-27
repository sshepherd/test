function elevation, ioff, phidiff, maxbeam, bmsep, bmnum, tfreq, tdiff, phi0, $
										debug=debug

C = 2.997e8		; speed of light in free space

; ioff[3] : interferometer offset [m]
;           +x along array towards increasing beam numbers
;           +y boresight direction with sign determined by x and z
;           +z up

; array separation in 3D
a = sqrt(ioff[1]*ioff[1] + ioff[0]+ioff[0] + ioff[2]*ioff[2])

; elevation angle correction, if arrays are at different heights [rad]
; note: that phidiff is the phase sign in hardware.dat file
;       it is +1 if cables are connected correctly, -1 if flipped
elev_corr = phidiff * asin(ioff[2]/a)

if (ioff[1] gt 0.) then begin	; interferometer in front of main
	phi_sign  = 1.
endif else begin							; interferometer behind main
	phi_sign  = -1.
  elev_corr = -elev_corr
endelse

; compute the angle (phi) the beam direction makes with the array boresite
offset = maxbeam/2. - 0.5
phi    = bmsep * (bmnum - offset) * !pi/180.
c_phi  = cos(phi)

; wave number in free space
k = 2*!pi * tfreq * 1000./C

; the phase difference phi0 is between -pi and +pi and gets positive,
; if the signal from the interferometer antenna arrives earlier at the
; receiver than the signal from the main antenna.
; If the cable to the interferometer is shorter than the one to
; the main antenna, than the signal from the interferometer
; antenna arrives earlier. tdiff < 0  --> dchi_cable > 0

dchi_path = - 2.*!pi * tfreq * 1000.0 * tdiff * 1.e-6

; If the interferometer antenna is in front of the main antenna
; then lower elevation angles correspond to earlier arrival
; and greater phase difference.
; If the interferometer antenna is behind of the main antenna
; then lower elevation angles correspond to later arrival
; and smaller phase difference

chi_max = phi_sign * k * a * c_phi + dchi_path

; change phi0 by multiples of twopi, until it is in the range
; (chi_max - twopi) to chi_max (interferometer in front)
; or chi_max to (chi_max + twopi) (interferometer in the back)

phi_temp = phi0 + 2.*!pi * floor( (chi_max - phi0)/(2.*!pi) )
if (phi_sign lt 0.) then phi_temp += 2.*!pi

; subtract the cable effect
psi = phi_temp - dchi_path

theta = psi/(k * a)
theta = (c_phi*c_phi - theta*theta)

; set elevation angle to 0 for out of range values

if ( (theta lt 0.) or (abs(theta) gt 1.) ) then begin
	if keyword_set(debug) then print, 'theta out of range: ',theta
	theta = -elev_corr	; ensure 0 at end
endif else begin
	theta = asin(sqrt(theta))
endelse

theta = 180. * (theta + elev_corr) / !pi	; not sure about the correction...
return, theta

;  double k;          /* wave number; 1/m */
;  double phi;        /* beam direction off boresight; rad */
;  double c_phi;      /* cosine of phi                     */
;  double dchi_cable; /* phase shift caused by cables; rad */
;  double chi_max;    /* maximum phase shift possible; rad */
;  double phi_temp;   /* actual phase angle + cable;   rad */
;  double psi;        /* actual phase angle - cable;   rad */
;  double theta;      /* angle of arrival for horizontal antennas; rad */
;  double offset=7.5; /* offset in beam widths to the edge of the array */
;  static double antenna_separation= 0.0; /* m */

end

