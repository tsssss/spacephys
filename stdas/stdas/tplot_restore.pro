;+
;PROCEDURE:  tplot_restore ,filenames=filenames, all=all, sort=sort
;PURPOSE:
;   Restores tplot data, limits, name handles, options, and settings.
;INPUT:
;KEYWORDS:
;   filenames:  file name or array of filenames to restore.  If
;               no file name is chosen and the all keyword is not set,
;		tplot_restore will look for and restore a file called
;		saved.tplot.
;   all: restore all *.tplot files in current directory
;   append: append saved data to existing tplot variables
;   sort: sort data by time after loading in
;   get_tvars: load tplot_vars structure (the structure containing tplot
;		options and settings even if such a structure already exists
;		in the current session.  The default is to only load these
;		if no such structure currently exists in the session.
;   restored_varnames=the tplot variable names for the restored data
;SEE ALSO:      "TPLOT_SAVE","STORE_DATA", "GET_DATA", "TPLOT"
;
;CREATED BY:    Peter Schroeder
;LAST MODIFICATION:     Added restore_varnames, 19-jun-2007, jmm
;                       Changed the obsolete IDL routine str_sep to
;                         strsplit. Also added additional output
;						  text.                   21-may-2008, cg
;                       Removed additional output text - Use dprint,debug=3  to restore text.   Nov 2008
;
;-
pro tplot_restore,filenames=filenames,all=all,append=append,sort=sort,$
	get_tvars=get_tvars,verbose=verbose, restored_varnames=restored_varnames

@tplot_com.pro

tplot_quant__define

idx = where(strpos(filenames,'.tplot') eq -1, cnt)
if cnt ne 0 then filenames[idx]+= '.tplot'

if keyword_set(all) then filenames = file_search('*.tplot')
if size(/type,filenames) ne 7 then $
  filenames = 'saved.tplot'
n = n_elements(filenames)
restored_varnames = ''
for i=0,n-1 do begin
  restore,filenames(i),/relaxed
  if keyword_set(tv) then begin
	  chkverb = where(tag_names(tv.options) eq 'VERBOSE',verbosethere)
	  if not verbosethere then begin
	  	optstruct = tv.options
  		setstruct = tv.settings
  		newopt = create_struct(optstruct,'VERBOSE',0)
  		tv = 0
  		tv = {options: newopt, settings: setstruct}
  		optstruct = 0
  		setstruct = 0
  		newopt = 0
  	endif
  endif
  if (n_elements(tplot_vars) eq 0) or keyword_set(get_tvars) then $
  	if keyword_set(tv) then tplot_vars = tv
  if keyword_set(dq) then begin
  	for j=0,n_elements(dq.name)-1 do begin
  		thisdq = dq(j)
  		dprint,dlevel=3, 'The tplot variable '+thisdq.name+' is being restored.'
        restored_varnames = [restored_varnames, thisdq.name]
  		names = strsplit(thisdq.name,'.')
  		if keyword_set(append) then get_data,thisdq.name,ptr=olddata
  		if keyword_set(append) and keyword_set(olddata) then begin
  		   if keyword_set(*thisdq.dh) then begin
  				if thisdq.dtype eq 1 then begin
 					newx = ptr_new([*olddata.x,*(*thisdq.dh).x])
					newy = ptr_new([*olddata.y,*(*thisdq.dh).y])
  					ptr_free,(*thisdq.dh).x,(*thisdq.dh).y
  					oldv = ptr_new()
  					str_element,olddata,'v',oldv
  					if ptr_valid(oldv) then begin
  						if ndimen(*oldv) eq 1 then $
  							newv = ptr_new(*oldv) else $
  							newv = ptr_new([*oldv,*(*thisdq.dh).v])
  						ptr_free,(*thisdq.dh).v
  						newdata={x: newx, y: newy, v: newv}
  					endif else newdata={x: newx, y: newy}
  					olddata = 0
  				endif else begin
  					newdata = olddata
  					dattags = tag_names(olddata)
  					for k = 0,n_elements(dattags)-1 do begin
  						str_element,*thisdq.dh,dattags(k),foo
  						foo = *foo
  						str_element,newdata,dattags(k),[*olddata.(k),foo],/add
  					endfor
  				endelse
  				store_data,verbose=verbose,thisdq.name,data=newdata
  		   endif
  		endif else begin
  		   store_data,verbose=verbose,thisdq.name,data=*thisdq.dh,limit=*thisdq.lh,dlimit=*thisdq.dl,/nostrsw
         dprint,dlevel=3, 'The tplot variable '+thisdq.name+' has been restored.'
  		endelse
  		if keyword_set(sort) then tplot_sort,thisdq.name
  	endfor
	ptr_free,dq.dh,dq.dl,dq.lh
  endif
  dq = 0
  tv = 0
endfor
If(n_elements(restored_varnames) Gt 1) Then $
  restored_varnames = restored_varnames[1:*]
end
