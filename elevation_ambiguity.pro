pro elevation_ambiguity

; 2 pi ambiguity plot

bmsep   = 3.24
maxbeam = 24
offset  = maxbeam/2. - 0.5

C = 2.997e8
phi = -offset*bmsep + findgen(maxbeam)*bmsep

y0 = 40.
yN = 150.
dy = 5.
ny = fix((yN-y0)/dy)+1

set_plot, 'PS'
device, file='2pi.ps', /color
loadct, 1, /silent

!x.margin=[6,9]
a=findgen(17)*!pi/8.
usersym, cos(a), sin(a), /fill
sz = .5

for m=0,ny-1 do begin
	; special case of Y offset only
	yoff    = y0 + m*dy
	ioff    = [0., -yoff, 0]		; offset
	sep     = abs(ioff[1])
	nfreq   = 11
	tfreq   = 8000. + findgen(nfreq)*1000		; 8-18 MHz

	for k=0,nfreq-1 do begin
		freq = tfreq[k]
		theta0 = acos((sep - C/(freq*1000.))/(cos(phi*!dtor)*sep))
		if k eq 0 then begin
			plot, phi, theta0/!dtor, /nodata, $
						xtit='Beam Direction !7u!3!L0!N (degrees)', $
						ytit='Elevation Angle !7h!3 (degrees)', $
						tit='Maximum Elevation Angle for 2!7p!3 Ambiguity', $
						yrange=[0,90], /ystyle
			xyouts, 0,20, strtrim(fix(yoff),2)+' m', align=.5
		endif

		oplot, phi, theta0/!dtor
		oplot, phi, theta0/!dtor, color=k*20, psym=8, symsize=sz
		xyouts, !x.crange[1], theta0[maxbeam-1]/!dtor, $
						' '+string(format='(f4.1)', freq/1000.)+ ' MHz'
	endfor
endfor
device, /close

end
