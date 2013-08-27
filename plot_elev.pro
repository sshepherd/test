; Testing for CVR elevation angles August 2103
;
; plot elevation angles v range for a single period
; reads data file e.out that was dumped using fanplot.pro with /dumpelev

pro plot_elev

id = 0 & rng = 0 & gflg = 0
bm = 0 & xqf = 0 & xgfl = 0
nr = 0

k = 1000
arng  = intarr(k,2)
aelv  = fltarr(k,2)
aphi0 = fltarr(k,2)
axqf  = intarr(k,2)
abm   = intarr(k,2)
nelv  = intarr(2)

file = 'e.out'
openr, fp, file, /get_lun
while not eof(fp) do begin
	readf, fp, id,bm,nr
	for k=0,nr-1 do begin
		readf, fp, rng, xqf, gflg, xgfl, phi0,phi0_e, elv,elv_l,elv_h
		arng[nelv[id],id]  = rng
		aelv[nelv[id],id]  = elv
		aphi0[nelv[id],id] = phi0
		axqf[nelv[id],id]  = xqf
		abm[nelv[id],id]   = bm
		nelv[id]++
	endfor
endwhile
free_lun, fp
print, nelv

; plotting
set_plot, 'PS'
device, file='cv_elev.ps', /color
device, ysize=27, yoff=1, xsize=21, xoff=1
!p.multi=[0,1,4]
!x.charsize=1.5
!y.charsize=1.5
;!y.margin=[0,.5]
;!y.omargin=[8,0]
!y.omargin[0]=4

loadct, 39, /silent
a=findgen(17)*!pi/8.
sz = .5

;el = 0.
;eh = 60.
;de = 0.5
;ne = fix((eh-el)/de) + 1

nrng = 110
maxe = 50
sites = ['CVE','CVW']
cols  = [50,254]
dcol  = 10
for k=0,1 do begin
;	if k eq 0 then begin
;		xtn = replicate(' ',30)
;		xt  = ''
;	endif else begin
		xtn = replicate('',30)
		xt  = 'Range Gate'
;	endelse
	plot, [0,70],[0,maxe], /nodata, xtickn=xtn, xtit=xt, $
				ytit='Elevation Angle (degrees)'
	xyouts, !x.crange[0],!y.crange[1], '!C!C  '+sites[k], charthick=2

	rng = arng[0:nelv[k]-1,k]
	elv = aelv[0:nelv[k]-1,k]
	xqf = axqf[0:nelv[k]-1,k]
	bm  = abm[0:nelv[k]-1,k]

	; good xcf
	q = where(xqf eq 1, nq)
	usersym, cos(a), sin(a), /fill
;	plots, rng[q], elv[q], psym=8, symsize=sz, color=cols[k]
	for m=0,nq-1 do $
		plots, rng[q[m]], elv[q[m]], psym=8, symsize=sz, color=bm[q[m]]*dcol

	; bad xcf
	q = where(xqf ne 1, nq)
	badx = intarr(nrng+1,5)
	usersym, cos(a), sin(a)
	for m=0,nq-1 do badx[rng[q[m]],xqf[q[m]]]++
	for m=0,nrng do $
		for n=0,4 do $
			if badx[m,n] gt 0 then $
				plots, m,n, psym=8, symsize=.1*badx[m,n], color=0
;				xyouts, m,n-.5, strtrim(badx[m,n],2), align=.5, charsize=.5
;	plots, rng[q], xqf[q], psym=8, symsize=sz*.75, color=0
;print, badx
endfor

for k=0,1 do begin
	if k eq 0 then begin
		xtn = replicate(' ',30)
		xt  = ''
	endif else begin
		xtn = replicate('',30)
		xt  = 'Range Gate'
	endelse
	plot, [0,70],[-1,1]*!pi, /nodata, xtickn=xtn, xtit=xt, $
				ytit='phi0 (radians)'
	xyouts, !x.crange[0],!y.crange[1], '!C!C  '+sites[k], charthick=2

	rng  = arng[0:nelv[k]-1,k]
	phi0 = aphi0[0:nelv[k]-1,k]
	xqf  = axqf[0:nelv[k]-1,k]
	bm   = abm[0:nelv[k]-1,k]

	; good xcf
	q = where(xqf eq 1, nq)
	usersym, cos(a), sin(a), /fill
;	plots, rng[q], phi0[q], psym=8, symsize=sz, color=cols[k]
	for m=0,nq-1 do $
		plots, rng[q[m]], phi0[q[m]], psym=8, symsize=sz, color=bm[q[m]]*dcol

	; bad xcf
	q = where(xqf ne 1, nq)
	badx = intarr(nrng+1,5)
	usersym, cos(a), sin(a)
	for m=0,nq-1 do badx[rng[q[m]],xqf[q[m]]]++
	for m=0,nrng do $
		for n=0,4 do $
			if badx[m,n] gt 0 then $
				plots, m,n, psym=8, symsize=.1*badx[m,n], color=0
endfor

device, /close
end
