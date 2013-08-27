pro test_elevation

; testing existing algorithms; elevation.c and elev_goose.c

ioff    = [0., -80., 0]		; CV offsets
maxbeam = 24
bmnum   = 12			; beam number
tfreq   = 11000.	; 11 MHz
phidiff = 1				; phase sign
tdiff   = -0.351	; in us
;tdiff   = 0.
bmsep   = 3.24

; computing what phi0 aught to be for this elevation angle
ang = 30. ; want this elevation angle
ang = 15. ; want this elevation angle
sep = sqrt(ioff[0]*ioff[0] + ioff[1]*ioff[1] + ioff[2]*ioff[2])
len = abs(ioff[1])*cos(ang*!dtor)
;print, len

C = 2.997e8
lambda = C / (tfreq * 1000.)
k = 2*!pi * tfreq * 1000./C
phi0 = fix(len*k/(2*!pi))*2*!pi
phi0 = phi0 - len*k ; ignoring tdiff here...
;print, len, lambda, k, phi0, phi0/(2*!pi)

theta = elevation(ioff,phidiff,maxbeam,bmsep,bmnum,tfreq,tdiff,phi0,/debug)
print, theta
;      29.9471 for tdiff = 0 and angle = 30
;      14.9082 for tdiff = 0 and angle = 15

theta = elev_goose(ioff,phidiff,maxbeam,bmsep,bmnum,tfreq,tdiff,phi0,/debug)
print, theta
;      infinite loop for tdiff = 0 ....

end
