;+
;PROCEDURE:	get_en_spec
;PURPOSE:	
;	Generates energy-time spectrogram data structures for tplot
;INPUT:		
;	data_str, 	a string (either 'eh','pl','eesa_surv','ess', ...)
;			where get_'string' returns a 2D or 3D 
;			data structure
;KEYWORDS:
;	T1:		start time, seconds since 1970
;	T2:		end time, seconds since 1970		
;	ANGLE:		fltarr(2),fltarr(4)	angle range to sum over
;	ARANGE:		intarr(2)		bin range to sum over
;	BINS:		bytarr(dat.nbins)	bins to sum over
;	gap_time: 	time gap big enough to signify a data gap 
;			(default 200 sec, 8 sec for FAST)
;	NO_DATA: 	returns 1 if no_data else returns 0
;	UNITS:		convert to these units if included
;	NAME:  		New name of the Data Quantity
;	BKG:  		A 3d data structure containing the background counts.
;	FLOOR:  	Sets the minimum value of any data point to sqrt(bkg).
;	MISSING: 	value for bad data.
;	RETRACE: 	Set to number of retrace energy steps to be eliminated starting at energy step 0
;	CALIB:		Calib keyword passed on to get_"get_dat"_ts
;
;
;
;CREATED BY:	J.McFadden
;VERSION:	1
;LAST MODIFICATION:  97/03/04
;MOD HISTORY:
;		97/03/04	T1,T2 keywords added
;		97/05/22	CALIB keyword added
;
;
;NOTES:	  
;	Current version only works for FAST
;-

pro get_en_spec,data_str,  $
;       ram_array, $
	T1=t1, $
	T2=t2, $
;	ENERGY=en, $
;	ERANGE=er, $
;	EBINS=ebins, $
	ANGLE=an, $
	ARANGE=ar, $
	BINS=bins, $
	gap_time=gap_time, $ 
	no_data=no_data, $
	units = units,  $
        name  = name, $
	bkg = bkg, $
        missing = missing, $
        floor = floor, $
        retrace = retrace, $
        CALIB = calib,  $
        no_ram = no_ram



;	Time how long the routine takes
ex_start=systime(1)

;	Set defaults for keywords, etc.

n = 0
max = 30000        ; this could be improved
all_same = 1

routine = 'get_'+data_str

;t = 0             ; get first sample
;dat = call_function(routine,t,CALIB=calib,index=0)
; FAST uses "start" keyword

if keyword_set(t1) then begin
	t=t1
	if routine eq 'get_fa_sebs' then dat = call_function(routine,t,/first) else dat = call_function(routine,t,CALIB=calib)
endif else begin
	t = 1000             ; get first sample
	dat = call_function(routine,t,CALIB=calib, /start)
endelse

if dat.valid eq 0 then begin no_data = 1 & return & end $
else no_data = 0

ytitle = data_str + '_en_spec'
last_time = (dat.time+dat.end_time)/2.
nbins = dat.nbins
nmaxvar = dat.nenergy

default_gap_time = 200.
if dat.project_name eq 'FAST' then begin
	nmaxvar=96
	default_gap_time = 8.
endif
if not keyword_set(gap_time) then gap_time = default_gap_time

time   = dblarr(max)
data   = fltarr(max,nmaxvar)
var   = fltarr(max,nmaxvar)
nvar = dat.nenergy
nmax=nvar

if not keyword_set(units) then units = 'Counts'
if not keyword_set(missing) then missing = !values.f_nan

;	Collect the data - Main Loop
    
; May want to use the following lines when "index" is operational in FAST get* routines    
;times = call_function(routine,t,CALIB=calib, /times)
;for idx=0,n_elements(times)-1 do begin
;    if (dat.valid eq 0) or (n ge max) then  goto, continue  ; goto YUCK!

if keyword_set(t2) then tmax=t2 else tmax=1.e30

while (dat.valid ne 0) and (n lt max) do begin
if (dat.valid eq 1) then begin

	count = dat.nbins
	if keyword_set(an) then bins=angle_to_bins(dat,an)
	if keyword_set(ar) then begin
		nb=dat.nbins
		bins=bytarr(nb)
		if ar(0) gt ar(1) then begin
			bins(ar(0):nb-1)=1
			bins(0:ar(1))=1
		endif else begin
			bins(ar(0):ar(1))=1
		endelse
	endif
; Set the "count" to the number of bins summed over
	if not keyword_set(bins) then ind=indgen(dat.nbins) else ind=where(bins,count)
	if units eq 'Counts' then norm = 1 else norm = count

	if abs((dat.time+dat.end_time)/2.-last_time) ge gap_time then begin
		if n ge 2 then dbadtime = time(n-1) - time(n-2) else dbadtime = gap_time/2.
		time(n) = (last_time) + dbadtime
		data(n,*) = missing
		var(n,*) = missing
		n=n+1
		if (dat.time+dat.end_time)/2. gt time(n-1) + gap_time then begin
			time(n) = (dat.time+dat.end_time)/2. - dbadtime
			data(n,*) = missing
			var(n,*) = missing
			n=n+1
		endif
	endif

	if keyword_set(bkg) then dat = sub3d(dat,bkg)
	if keyword_set(units) then dat = conv_units(dat,units)

	nvar = dat.nenergy
	if nvar gt nmax then nmax = nvar
	time(n)   = (dat.time+dat.end_time)/2.
	if ind(0) ne -1 then begin
		data(n,0:nvar-1) = total( dat.data(*,ind), 2)/norm
		var(n,0:nvar-1) = total( dat.energy(*,ind), 2)/count
	endif else begin
		data(n,0:nvar-1) = 0
		var(n,0:nvar-1) = total( dat.energy(*,0), 2)
	endelse

; test the following lines, the 96-6-19 version of tplot did not work with !values.f_nan
;	if nvar lt nmaxvar then data(n,nvar:nmaxvar-1) = !values.f_nan
;	if nvar lt nmaxvar then var(n,nvar:nmaxvar-1) = !values.f_nan
	if nvar lt nmaxvar then data(n,nvar:nmaxvar-1) = data(n,nvar-1)
	if nvar lt nmaxvar then var(n,nvar:nmaxvar-1) = 1.5*var(n,nvar-1)-.5*var(n,nvar-2)

	if (all_same eq 1) then begin
		if dimen1(where(var(n,0:nvar-1) ne var(0,0:nvar-1))) gt 1 then all_same = 0
	endif
	last_time = time(n)
	n=n+1

endif else begin
	print,'Invalid packet, dat.valid ne 1, at: ',time_to_str(dat.time)
endelse

	dat = call_function(routine,t,CALIB=calib,/ad)
;	dat = call_function(routine,t,CALIB=calib,index=idx)
	if dat.valid ne 0 then if dat.time gt tmax then dat.valid=0

endwhile
;endfor
;continue:

;	Store the data

	if count ne nbins then ytitle = ytitle+'_'+strtrim(count,2)
	if keyword_set(name) eq 0 then name=ytitle else ytitle = name
	ytitle = ytitle+' ('+units+')'

if not keyword_set(retrace) then begin
;	If you want to plot the retrace, set the retrace flag to 1.
	data = data(0:n-1,0:nmax-1)
	var = var(0:n-1,0:nmax-1)
endif else begin
	data = data(0:n-1,retrace:nmax-1)
	var = var(0:n-1,retrace:nmax-1)
endelse

print,'all_same=',all_same
;labels=''
; The following has be removed so that FAST summary cdf files contain both arrays
;if all_same then begin
;	var = reform(var(0,*))
;	labels = strtrim( round(var) ,2)+ ' eV'
;endif

time = time(0:n-1)

if keyword_set(t1) then begin
	ind=where(time ge t1)
	time=time(ind)
	data=data(ind,*)
	var=var(ind,*)
endif
if keyword_set(t2) then begin
	ind=where(time le t2)
	time=time(ind)
	data=data(ind,*)
	var=var(ind,*)
endif

;datastr = {ztitle:units,x:time,y:data,v:var,  $
;	labels:labels,	$
;    ylog:1,panel_size:2.}

;; if keyword_set(no_ram) then data=remove_ram(data,ram_array=ram_array)
if keyword_set(no_ram) then data=remove_ram(data)

datastr = {x:time,y:data,v:var}
store_data,name,data=datastr
ex_time = systime(1) - ex_start
message,string(ex_time)+' seconds execution time.',/cont,/info
print,'Number of data points = ',n

return

end
