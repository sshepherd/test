pro ebars, x,y, dx=dx, dy=dy
	if not keyword_set(dx) then dx = 1
	if not keyword_set(dy) then dy = 1

	nx = n_elements(x)
	for k=0,nx-1 do begin
		plots, x[k]+[0,0], y[k]+[-1,1]*dy
		plots, x[k]+[-1,1]*dx, y[k]+[1,1]*dy
		plots, x[k]+[-1,1]*dx, y[k]-[1,1]*dy
	endfor
end

pro plot_2D_data, files

nf   = n_elements(files)
k    = 30
alph = fltarr(nf,k)
azm  = fltarr(nf,k)
gain = fltarr(nf,k)
bw   = fltarr(nf,k)
nval = intarr(nf)
tfrq = fltarr(nf)
aaz0 = fltarr(nf)
bet  = fltarr(nf,k)

for m=0,nf-1 do begin
	openr, fp, files[m], /get_lun
	readf, fp, freq, az0
	k = 0
	while not eof(fp) do begin
		readf, fp, el, az, gn, w
		alph[m,k] = el
		azm[m,k]  = az
		gain[m,k] = gn
		bw[m,k]   = w
		k++
	endwhile
	free_lun, fp

	nval[m] = k
	tfrq[m] = freq
	aaz0[m] = az0

	; Ray's formula
	calp  = cos(alph[m,*]*!dtor)	; cosine of elevation angle
	cbet0 = cos((90-az0)*!dtor)		; cos of comp. at 0 elev angle (beam angle)
	bet[m,*] = 90- acos(cbet0/calp)/!dtor

endfor

set_plot, 'PS'
device, file='elevation_data.ps', /color
loadct, 39, /silent

a = findgen(17)*!pi/8.
usersym, cos(a), sin(a), /fill

plot, alph[0,*], azm[0,*], /nodata, $
			xrange=[0,50], xtit='Elevation Angle (degrees)', $
			yrange=[10,45], ytit='Beam Azimuth (degrees)'

for m=0,nf-1 do begin
	oplot, !x.crange, aaz0[m]+[0,0], lines=5

	oplot, alph[m,*],azm[m,*], psym=8, color=254-m*50, symsize=1-.2*m
	ebars, alph[m,*],azm[m,*], dx=.25, dy=.5

	oplot, alph[m,0:nval[m]-1], bet[m,0:nval[m]-1], /lines

;xyouts, .5*total(!x.crange), !y.crange[1], $
;				'!C!C'+string(format='(f4.1)', freq*1e-3)+' MHz', align=.5

endfor
device, /close
end
