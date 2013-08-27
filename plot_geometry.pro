;-----------------------------------------------------------------------------
;
; show examples of new way using distance between two planes, comparing
;  corrected, uncorrected and current algorithms
;

pro plot_geometry, psname=psname, tfreq=tfreq, phidiff=phidiff, tdiff=tdiff, $
										zparm=zparm, phi_obs=phi_obs, beams=beams

if not keyword_set(psname) then psname='new.ps'
set_plot, 'PS'
device, file=psname, /color
device, ysize=26, yoff=1, xsize=21, xoff=1
loadct, 39, /silent
!x.margin[1] = 6

if not keyword_set(tfreq) then   tfreq   = 10000.	; TX freq in kHz
if not keyword_set(phidiff) then phidiff = 1			; phase sign
if not keyword_set(tdiff) then   tdiff   = 0.

bmsep   = 3.24
maxbeam = 24
if not keyword_set(beams) then beams = [7,8,9,10,11,12,13,14,15,16]

Y = -80.

x0 = -50
xN =  50
dx = 1.
nx = fix((xN-x0)/dx)+1

if not keyword_set(zparm) then begin
	X = x0 + findgen(nx)*dx		; vary X
	Z = 0.	; Z is fixed
	xz = X
	xtit = 'X offset (m)'
endif else begin
	Z = x0 + findgen(nx)*dx		; vary Z
	X = 0.	; X is fixed
	xz = Z
	xtit = 'Z offset (m)'
endelse

offset  = maxbeam/2. - 0.5
phi0s   = bmsep * (beams - offset); * !pi/180.
nphi  = n_elements(phi0s)

; range of observed delays
if keyword_set(phi_obs) then phi_obs_saf = phi_obs else phi_obs_saf = .8*!pi

cz = .75
plot, xz,xz, /nodata, $
				xtit=xtit, /xstyle, $
				ytit='Elevation Angle [degrees]', yrange=[0,70]

plots, [0,0],!y.crange, /lines
theta = fltarr(nx,3)

for ii=0,nphi-1 do begin
	bmnum   = beams[ii]

	for i=0,nx-1 do begin

		if not keyword_set(zparm) then ioff = [X[i],Y,Z] else ioff = [X,Y,Z[i]]

phi_obs = phi_obs_saf
		theta[i,0] = elev_new(ioff,phidiff,maxbeam,bmsep,bmnum,tfreq,tdiff,phi_obs)
phi_obs = phi_obs_saf
		theta[i,1] = elev_new(ioff,phidiff,maxbeam,bmsep,bmnum,tfreq,tdiff,phi_obs,$
									/corr)
phi_obs = phi_obs_saf
		theta[i,2] = elevation(ioff,phidiff,maxbeam,bmsep,bmnum,tfreq,tdiff,phi_obs)
	endfor

;	plots, X,theta[*,0], color=75, thick=2		; new, no correction
	plots, xz,theta[*,1], color=50, thick=2		; new, correction
	if phi0s[ii] lt 0 then begin
		xyouts, !x.crange[0], theta[0,1], string(format='(i2.2)', bmnum), $
						/align, charsize=cz, color=50
	endif else begin
		xyouts, !x.crange[1], theta[nx-1,1], string(format='(i2.2)', bmnum), $
						charsize=cz, color=50
	endelse
	plots, xz,theta[*,2], color=250, thick=4, /lines	; current
	if phi0s[ii] lt 0 then begin
		xyouts, !x.crange[0], theta[0,2], string(format='(i2.2)', bmnum)+'   ', $
						/align, charsize=cz, color=250
	endif else begin
		xyouts, !x.crange[1], theta[nx-1,2], '   '+string(format='(i2.2)', bmnum), $
						charsize=cz, color=250
	endelse

	xyouts, !x.crange[0], !y.crange[1], '!C!C!C  tfreq '+ $
						string(format='(f4.1)', tfreq*1e-3)+' MHz' + $
						'!C  delay '+ string(format='(f5.2)', phi_obs) + ' radians'

	if not keyword_set(zparm) then begin
		xyouts, !x.crange[1], !y.crange[1], '!C!C!C[Y,Z] = ['+$
					string(format='(f5.1,",",f5.1)', ioff[1:2])+']  !C'+$
					'tdiff = '+string(format='(i4)', round(tdiff*1e3))+' ns  ', /align
	endif else begin
		xyouts, !x.crange[1], !y.crange[1], '!C!C!C[X,Y] = ['+$
					string(format='(f5.1,",",f5.1)', ioff[0:1])+']  !C'+$
					'tdiff = '+string(format='(i4)', round(tdiff*1e3))+' ns  ', /align
	endelse

;	; plot differences
;	plot, delays, delays, /ylog, yrange=[1e-8,1e2], /nodata, /ystyle, $
;				xtit='Delay !7h!3 (radians)', $
;				ytit='Difference'
;
;	plots, delays, abs(theta[*,1]-theta[*,2]), color=50
;	plots, delays, abs(theta[*,0]-theta[*,2]), color=25

endfor

device, /close

end

