;-----------------------------------------------------------------------------
;
; show examples of new way using distance between two planes, comparing
;  corrected, uncorrected and current algorithms
;

pro simulate_elev, psname=psname, ioff=ioff, maxbeam=maxbeam, bmsep=bmsep, $
									tfreq=tfreq, phidiff=phidiff, tdiff=tdiff, beams=beams, $
									yr=yr

if not keyword_set(psname) then psname='new.ps'
set_plot, 'PS'
device, file=psname, /color
device, ysize=26, yoff=1, xsize=21, xoff=1
loadct, 39, /silent
!p.multi=[0,1,2]
!x.margin[1]=12
!y.margin=[0,.5]
!y.omargin=[6,.5]

if not keyword_set(ioff) then    ioff    = [0., -90., 0]	; 3 lambda @ 10 MHz
if not keyword_set(maxbeam) then maxbeam = 24
if not keyword_set(bmsep) then   bmsep   = 3.24
if not keyword_set(tfreq) then   tfreq   = 10000.	; TX freq in kHz
if not keyword_set(phidiff) then phidiff = 1			; phase sign
if not keyword_set(tdiff) then   tdiff   = 0.
if not keyword_set(beams) then   beams   = maxbeam/2.-.5	; boresight
if not keyword_set(yr) then yr = [-20,0]

X = ioff[0]
Y = ioff[1]
Z = ioff[2]

offset  = maxbeam/2. - 0.5
phi0s   = bmsep * (beams - offset); * !pi/180.
nphi  = n_elements(phi0s)

print, phi0s
tlow  =   0		; range of thetas (elevation angles)
thgh  =  80
tdel  = .1
nthe  =  fix((thgh-tlow)/tdel)+1
theta = tlow + findgen(nthe)*tdel

plot, theta, theta, /nodata, $
			xtit='Elevation Angle !7h!3 (degrees)', $
			ytit='Delay !7u!L!3total!N [radians]', yrange=yr

plots, !x.crange,[0,0], /lines

xyouts, !x.crange[0],!y.crange[1], 'beam   ', /align, charthick=2
for i=0,nphi-1 do begin
	;+--------------------------------------------------------------------------
	; beam angle in horizontal plane (x-y)
	;+--------------------------------------------------------------------------

	phi0 = phi0s[i]
	print, 'beam: ',beams[i]
	print, 'phi0: ',phi0

	; no correction
	phi = phi0
; actual phi (beam direction) using Ray's formula
	phi_a = asin(sin(phi0*!dtor)/cos(theta*!dtor))/!dtor

	;+--------------------------------------------------------------------------
	; distance between two planes
	;+--------------------------------------------------------------------------

	; no correction
	d_old = (X*sin(phi*!dtor) + Y*cos(phi*!dtor))*cos(theta*!dtor) + $
					Z*sin(theta*!dtor)
	; correction
	d = (X*sin(phi_a*!dtor) + Y*cos(phi_a*!dtor))*cos(theta*!dtor) + $
					Z*sin(theta*!dtor)

	C = 2.997e8
	k = 2*!pi * tfreq * 1000./C

	;+--------------------------------------------------------------------------
	; phase delay [radians] due to array geometry
	;+--------------------------------------------------------------------------

	phi_geo   = k*d					; correction
	phi_geo_old = k*d_old		; no correction

	;+--------------------------------------------------------------------------
	; phase delay [radians] due to electrical path difference
	;+--------------------------------------------------------------------------
	phi_path  = 2.*!pi * tfreq*1000. * tdiff*1e-6

	;+--------------------------------------------------------------------------
	; total phase delay [radians] due to electrical path and geometry
	;+--------------------------------------------------------------------------

	phi_total     = phi_geo - phi_path					; correction
	phi_total_old = phi_geo_old - phi_path			; no correction

;	print, 'phi0      = ', phi0
;	print, 'phi_geo   = ', phi_geo
;	print, 'phi_path  = ', phi_path
;	print, 'phi_total = ', phi_total

	;+--------------------------------------------------------------------------
	; making up a reasonable phi_obs; this is innaccurate
	;+--------------------------------------------------------------------------
	ang = 0 ; want this elevation angle
	ang = 50. ; want this elevation angle (above 2pi ambiguity!!)
	ang = 30. ; want this elevation angle; only approximate
	ang = 47. ; want this elevation angle; only approximate

	len      = -Y*cos(ang*!dtor)
	npi      = fix((len*k+phi_path)/(2*!pi))
	phi_obs  = len*k + phi_path	; estimate of phi_obs
;print, phi_obs/(2*!pi)
	n2pi     = fix(phi_obs/(2*!pi))
;print, n2pi
	phi_obs -= n2pi*2*!pi
;print, phi_obs/(2*!pi)
;stop
	phi_obs = -phi_obs

	;max_delay = max(phi_total)		; will be min if flipped?
	;d2pi = (fix((phi0-max_delay)/(2*!pi))+1)*2*!pi		; flip?
	n2pi = fix((phi_obs-phi_total[0])/(2*!pi))		; flip?
	n2pi_old = fix((phi_obs-phi_total_old[0])/(2*!pi))		; flip?
	d2pi = n2pi*2*!pi		; flip?
	d2pi_old = n2pi_old*2*!pi		; flip?

	;;print, 'phi_obs = ', phi_obs
	;print, 'max_delay= ', max_delay
	;;print, 'd2pi    = ', d2pi

	phi_obs_saf = phi_obs
	;+--------------------------------------------------------------------------
	; uncorrected solution
	;+--------------------------------------------------------------------------
	plots, theta, phi_total_old, color=75, thick=4
	xyouts, !x.crange[0], phi_total_old[0], $
					string(format='(i2.2)', beams[i])+'   ', /align, charthick=2
	oplot, !x.crange, phi_obs+[0,0], color=225, thick=2, /lines
	xyouts, !x.crange[1],phi_obs, ' original!C delay !7u!L!3obs.!N'
	phi_obs -= d2pi_old
	if n2pi_old lt 0 then sgn = '+' else sgn = '-'
	xyouts, !x.crange[1],phi_obs, $
					' !7u!L!3obs.!N'+sgn+strtrim(abs(2*n2pi_old),2)+'!7p!3!N'

	; now solve for theta
	D = X*sin(phi*!dtor) + Y*cos(phi*!dtor)
	E = (phi_path + phi_obs)/k
	sqrtfac = D*D*(Z*Z + D*D - E*E)
	denom = Z*Z + D*D
	nump  = E*Z + sqrt(sqrtfac)
	numn  = E*Z - sqrt(sqrtfac)

	tnocorr = asin(nump/denom)/!dtor
;	tnocorr = abs(asin(numn/denom)/!dtor)
	print, 'uncorrect: ',tnocorr
	oplot, tnocorr+[0,0], [!y.crange[0],phi_obs], /lines, color=225, thick=2
	oplot, [tnocorr,!x.crange[1]], phi_obs+[0,0], color=225, thick=2

; 2pi ambiguity, but specific for -Y only
;theta0 = acos((sep - C/(tfreq*1000.))/(cos(phi*!dtor)*sep))

	;+--------------------------------------------------------------------------
	; corrected solution
	;+--------------------------------------------------------------------------
	phi_obs = phi_obs_saf		; reset phi_obs
	plots, theta, phi_total, color=50, thick=4
	phi_obs -= d2pi
	if n2pi lt 0 then sgn = '+' else sgn = '-'
	xyouts, !x.crange[1],phi_obs, $
					' !7u!L!3obs.!N'+sgn+strtrim(abs(2*n2pi),2)+'!7p!3!N'

	; now solve for theta
	D = 1. - sin(phi0*!dtor)*sin(phi0*!dtor)
	E = (phi_path + phi_obs)/k - sin(phi0*!dtor)*X
	sqrtfac = -E*E + (Z*Z+Y*Y)*D
	denom = Z*Z + Y*Y
	nump  = E*Z + abs(Y)*sqrt(sqrtfac)
	numn  = E*Z - abs(Y)*sqrt(sqrtfac)

	tcorr = asin(nump/denom)/!dtor
;	tcorr = abs(asin(numn/denom)/!dtor)
	print, 'corrected: ',tcorr
	oplot, tcorr+[0,0], [!y.crange[0],phi_obs], /lines, color=254, thick=2
	oplot, [0,tcorr], phi_obs+[0,0], color=254, thick=2

	;+--------------------------------------------------------------------------
	; our current algorithem
	;+--------------------------------------------------------------------------
	phi_obs = phi_obs_saf		; reset phi_obs
	te = elevation(ioff,phidiff,maxbeam,bmsep,beams[i],tfreq,tdiff,phi_obs,/debug)
	oplot, te +[0,0], !y.crange, /lines, color=200, thick=2
	print, 'current:   ',te
endfor

;-----------------------------------------------------------------------------
;
; Show differences between techniques
;

; range of observed delays
d0 = -!pi 
dN = !pi
dd = !pi/101.
nd = fix((dN-d0)/dd) + 1
delays = d0 + findgen(nd)*dd

theta = fltarr(nd,3)

!p.multi=[0,1,2]
;!y.margin=[0,.5]
;!y.omargin=[6,0]

for ii=0,nphi-1 do begin
	bmnum   = beams[ii]

	for i=0,nd-1 do begin
		phi_obs = delays[i]
		theta[i,0] = elev_new(ioff,phidiff,maxbeam,bmsep,bmnum,tfreq,tdiff,phi_obs)
		theta[i,1] = elev_new(ioff,phidiff,maxbeam,bmsep,bmnum,tfreq,tdiff,phi_obs,$
									/corr)
		theta[i,2] = elevation(ioff,phidiff,maxbeam,bmsep,bmnum,tfreq,tdiff,phi_obs)
	endfor

	; plot elevation angles
	plot, delays, delays, /nodata, xtickn=replicate(' ',30), $
				ytit='Elevation Angle [degrees]', yrange=[0,60]

	plots, delays,theta[*,0], color=75, thick=2		; new, no correction
	plots, delays,theta[*,1], color=50, thick=2		; new, correction
	plots, delays,theta[*,2], color=200, thick=4, /lines	; current

	xyouts, !x.crange[0], !y.crange[1], '!C!C  beam '+$
						string(format='(i2.2)', bmnum)+'!C  tfreq '+$
						string(format='(f4.1)', tfreq*1e-3)+' MHz'

	xyouts, !x.crange[1], !y.crange[1], '!C!C[X,Y,Z] = ['+$
					string(format='(f5.1,",",f5.1,",",f5.1)', ioff)+']  !C'+$
					'tdiff = '+string(format='(i4)', round(tdiff*1e3))+' ns  ', /align

	; plot differences
	plot, delays, delays, /ylog, yrange=[1e-8,1e2], /nodata, /ystyle, $
				xtit='Delay !7h!3 (radians)', $
				ytit='Difference'

	plots, delays, abs(theta[*,1]-theta[*,2]), color=50
	plots, delays, abs(theta[*,0]-theta[*,2]), color=25

endfor

device, /close

end

