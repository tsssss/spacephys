;+
;PROCEDURE:   tplot  [,datanames]
;PURPOSE:
;   Creates a time series plot of user defined quantities.
;INPUT:
;   datanames: A string of space separated datanames.
;             wildcard expansion is supported.
;             if datanames is not supplied then the last values are used.
;             Each name should be associated with a data quantity.
;             (see the "STORE_DATA" and "GET_DATA" routines.)
;             Alternatively datanames can be an array of integers or strings.
;             run "TPLOT_NAMES" to show the current numbering.
;
;KEYWORDS:
;   TITLE:    A string to be used for the title. Remembered for future plots.
;   ADD_VAR:  Set this variable to add datanames to the previous plot.  If set
;         to 1, the new panels will appear at the top (position 1) of the
;         plot.  If set to 2, they will be inserted directly after the
;         first panel and so on.  Set this to a value greater than the
;         existing number of panels in your tplot window to add panels to
;             the bottom of the plot.
;   LASTVAR:  Set this variable to plot the previous variables plotted in a
;         TPLOT window.
;   PICK:     Set this keyword to choose new order of plot panels
;             using the mouse.
;   WINDOW:   Window to be used for all time plots.  If set to -1, then the
;             current window is used.
;   VAR_LABEL:  String [array]; Variable(s) used for putting labels along
;     the bottom. This allows quantities such as altitude to be labeled.
;   VERSION:  Must be 1,2,3, or 4 (3 is default)  Uses a different labeling
;   scheme.  Version 4 is for rocket-type time scales.
;   OVERPLOT: Will not erase the previous screen if set.
;   NAMES:    The names of the tplot variables that are plotted.
;   NOCOLOR:  Set this to produce plot without color.
;   TRANGE:   Time range for tplot.
;   NEW_TVARS:  Returns the tplot_vars structure for the plot created. Set
;         aside the structure so that it may be restored using the
;             OLD_TVARS keyword later. This structure includes information
;             about various TPLOT options and settings and can be used to
;             recreates a plot.
;   OLD_TVARS:  Use this to pass an existing tplot_vars structure to
;     override the one in the tplot_com common block.
;   GET_PLOT_POSITION: Returns an array containing the corners of each panel in the plot, to make it easier to overplot and annotate plots
;   HELP:     Set this to print the contents of the tplot_vars.options
;         (user-defined options) structure.
;
;RESTRICTIONS:
;   Some data must be loaded prior to trying to plot it.  Try running
;   "_GET_EXAMPLE_DAT" for a test.
;
;EXAMPLES:  (assumes "_GET_EXAMPLE_DAT" has been run)
;   tplot,'amp slp flx2' ;Plots the named quantities
;   tplot,'flx1',/ADD          ;Add the quantity 'flx1'.
;   tplot                      ;Re-plot the last variables.
;   tplot,var_label=['alt']   ;Put Distance labels at the bottom.
;       For a long list of examples see "_TPLOT_EXAMPLE"
;
;OTHER RELATED ROUTINES:
;   Examples of most usages of TPLOT and related routines are in
;      the crib sheet: "_TPLOT_EXAMPLE"
;   Use "TNAMES" function to return an array of current names.
;   Use "TPLOT_NAMES" to print a list of acceptable names to plot.
;   Use "TPLOT_OPTIONS" for setting various global options.
;   Plot limits can be set with the "YLIM" procedure.
;   Spectrogram limits can be set with the "ZLIM" procedure.
;   Time limits can be set with the "TLIMIT" procedure.
;   The "OPTIONS" procedure can be used to set all IDL
;      plotting keyword parameters (i.e. psym, color, linestyle, etc) as well
;      as some keywords that are specific to tplot (i.e. panel_size, labels,
;      etc.)  For example, to change the relative panel width for the quantity
;      'slp', run the following:
;            OPTIONS,'slp','panel_size',1.5
;   TPLOT calls the routine "SPECPLOT" to make spectrograms and
;      calls "MPLOT" to make the line plots. See these routines to determine
;      what other options are available.
;   Use "GET_DATA" to retrieve the data structure (or
;      limit structure) associated with a TPLOT quantity.
;   Use "STORE_DATA" to create new TPLOT quantities to plot.
;   The routine "DATA_CUT" can be used to extract interpolated data.
;   The routine "TSAMPLE" can also be used to extract data.
;   Time stamping is performed with the routine "TIME_STAMP".
;   Use "CTIME" or "GETTIME" to obtain time values.
;   tplot variables can be stored in files using "TPLOT_SAVE" and loaded
;      again using "TPLOT_RESTORE"
;
;CREATED BY:    Davin Larson  June 1995
;
;QUESTIONS?
;   See the archives at:  http://lists.ssl.berkeley.edu/mailman/listinfo/tplot
;Still have questions:
;   Send e-mail to:  tplot@ssl.berkeley.edu    someone might answer!
;
; $LastChangedBy: pcruce $
; $LastChangedDate: 2014-01-31 17:32:47 -0800 (Fri, 31 Jan 2014) $
; $LastChangedRevision: 14111 $
; $URL: svn+ssh://thmsvn@ambrosia.ssl.berkeley.edu/repos/spdsoft/trunk/general/tplot/tplot.pro $
;-

pro tplot,datanames,      $
   WINDOW = wind,         $
   NOCOLOR = nocolor,     $
   VERBOSE = verbose,     $
   wshow = wshow,         $
   OPLOT = oplot,         $
   OVERPLOT = overplot,   $
   TITLE = title,         $ ; sheng.
   LASTVAR = lastvar,     $
   ADD_VAR = add_var,     $
   LOCAL_TIME= local_time,$
   REFDATE = refdate,     $
   VAR_LABEL = var_label, $ ; sheng.
   OPTIONS = opts,        $
   T_OFFSET = t_offset,   $
   TRANGE = trng,         $
   NAMES = names,         $
   PICK = pick,           $
   new_tvars = new_tvars, $
   old_tvars = old_tvars, $
   datagap = datagap,     $  ;It looks like this keyword isn't actually used.  pcruce 10/4/2012
   get_plot_position=pos, $
   help = help, $
   novtitle = novtitle, $; sheng.
   noerase = noerase, $  ; sheng.
   uttick = uttick,$     ; sheng.
   nouttick = nouttick, $; sheng.
   short_vlabel = short_vlabel , $  ; sheng.
   figlabel = figlabs, $ ; sheng.
   xstyle = xstyle, $   ; sheng.
   ystyle = ystyle, $   ; sheng.
   vlab_margin = vlab_marg, $ ; sheng.
   position = pos0       ; sheng.


compile_opt idl2

@tplot_com.pro

; sheng. force no title, no label.
if n_elements(title) eq 0 then title = ''
if n_elements(var_label) eq 0 then var_label = ''

if 1 then begin   ; check for embedded calls
    stack = scope_traceback(/structure)
    stack = stack[0:n_elements(stack)-2]
    nocallsfrom = ['CTIME','TPLOT']
    incommon = array_union(nocallsfrom,stack.routine)
    w = where(incommon ne -1,nw)
    if nw gt 0 then begin
        dprint,dlevel=2,'Calls to TPLOT are not allowed from within '+strupcase(nocallsfrom[w])
        return
    endif
endif

if size(verbose,/type) eq 0 then  str_element,tplot_vars,'options.verbose',verbose ; get default verbose value if it exists
if size(wshow,/type) eq 0 then str_element,tplot_vars,'options.wshow',wshow


if keyword_set(old_tvars) then tplot_vars = old_tvars

if keyword_set(help) then begin
    printdat,tplot_vars.options,varname='tplot_vars.options'
    new_tvars = tplot_vars
    return
endif

; setup tplot_vars....
tplot_options,title=title,var_label=var_label,refdate=refdate, wind=wind, options = opts


if keyword_set(overplot) then oplot=overplot
if n_elements(trng) eq 2 then trange = time_double(trng)

chsize = !p.charsize
if chsize eq 0. then chsize=1.

def_opts= {ymargin:[4.,2.],xmargin:[12.,12.],position:fltarr(4), $
   title:'',ytitle:'',xtitle:'', $
   xrange:dblarr(2),xstyle:1,    $
   version:3, window:-1, wshow:0,  $
   charsize:chsize,noerase:0,overplot:0,spec:0}

extract_tags,def_opts,tplot_vars.options

; Define the variables to be plotted:

;str_element,tplot_vars,'options.varnames',tplot_var
; if n_elements(tplot_var) eq 0 then $
;    str_element,tplot_vars,'options.varnames',['NULL'],/add_replace

if keyword_set(pick) then $
   ctime,prompt='Click on desired panels. (button 3 to quit)',panel=mix,/silent
if n_elements(mix) ne 0 then datanames = tplot_vars.settings.varnames[mix]

if keyword_set(add_var)  then begin
   names = tnames(datanames,/all)
   if add_var eq 1 then datanames = [names,tplot_vars.options.varnames] else $
    if (add_var gt n_elements(tplot_vars.options.varnames)) then $
        datanames = [tplot_vars.options.varnames,names] else $
        datanames = [tplot_vars.options.varnames[0:add_var-2],names,$
           tplot_vars.options.varnames[add_var-1:*]]
endif


dt = size(/type,datanames)
ndim = size(/n_dimen,datanames)

if dt ne 0 then begin
   if dt ne 7 or ndim ge 1 then dnames = strjoin(tnames(datanames,/all),' ') $
   else dnames=datanames
endif else begin
	tpv_opt_tags = tag_names( tplot_vars.options)
	idx = where( tpv_opt_tags eq 'DATANAMES', icnt)
	if icnt gt 0 then begin
		dnames=tplot_vars.options.datanames
	endif else begin
		return
	endelse
endelse

;if dt ne 0 then names= tnames(datanames,/all)

if keyword_set(lastvar) then str_element,tplot_vars,'settings.last_varnames',names

;if keyword_set(names) then begin
;   str_element,tplot_vars,'settings.last_varnames',tplot_vars.options.varnames,$
;       /add_replace
;   str_element,tplot_vars,'options.varnames',names,/add_replace ;  array of names
;   str_element,tplot_vars,'settings.varnames',names,/add_replace
;endif else names = tplot_vars.options.varnames

str_element,tplot_vars,'options.lazy_ytitle',lazy_ytitle

varnames = tnames(dnames,nd,ind=ind,/all)

str_element,tplot_vars,'options.datanames',dnames,/add_replace
str_element,tplot_vars,'options.varnames',varnames,/add_replace

if nd eq 0 then begin
   dprint,dlevel=0,verbose=verbose,'No valid variable names found to tplot! (use TPLOT_NAMES to display)'
   return
endif

;ind = array_union(tplot_vars.options.varnames,data_quants.name)

sizes = fltarr(nd)
for i=0,nd-1 do begin
   dum = 1.
   lim = 0
   get_data,tplot_vars.options.varnames[i],alim=lim
   str_element,lim,'panel_size',value=dum
   sizes[i] = dum
endfor

plt = {x:!x,y:!y,z:!z,p:!p}

if (!d.flags and 256) ne 0  then begin    ; windowing devices
   current_window= !d.window > 0
   if def_opts.window ge 0 then w = def_opts.window $
   else w = current_window
;test to see if this window exists before wset, jmm, 7-may-2008:
;removed upper limit on window number, jmm, 19-mar-2009
   device, window_state = wins
   if(w Eq 0 Or wins[w]) then wset,w else begin
     dprint,verbose=verbose, 'Window is closed and Unavailable, Returning'
     w = current_window
     def_opts.window = w
     tplot_options, window = w
     return
   endelse
   if def_opts.wshow ne 0 || keyword_set(wshow) then wshow ;,icon=0   ; The icon=0 option doesn't work with windows
   str_element,def_opts,'wsize',value = wsize
   wi,w,wsize=wsize
endif

str_element,tplot_vars,'settings.y',replicate(!y,nd),/add_replace
str_element,tplot_vars,'settings.clip',lonarr(6,nd),/add_replace
str_element,def_opts,'ygap',value = ygap
str_element,def_opts,'charsize',value = chsize

if keyword_set(nocolor) then str_element,def_opts,'nocolor',nocolor,/add_replace

nvlabs = [0.,0.,0.,1.,0.,0.]
str_element,tplot_vars,'options.var_label',var_label
if keyword_set(var_label) then if size(/type,var_label) eq 7 then $
    if ndimen(var_label) eq 0 then var_label=tnames(var_label) ;,/extrac)
;ensure the index does not go out of range, other values will use default
if def_opts.version lt 1 or def_opts.version gt 5 then def_opts.version = 3
nvl = n_elements(var_label) + nvlabs[def_opts.version]
def_opts.ymargin = def_opts.ymargin + [nvl,0.]

!p.multi = 0
pos = keyword_set(pos0)? pos0: $    ; sheng.
    plot_positions(ysizes=sizes,options=def_opts,ygap=ygap)

if  keyword_set(trange) then str_element,tplot_vars,'options.trange',trange,/add_replace $
else  str_element,tplot_vars,'options.trange',trange
if trange[0] eq trange[1] then $
    trg=minmax(reform(data_quants[ind].trange),min_value=0.1) $
else trg = trange

tplot_var_labels,def_opts,trg,var_label,local_time,pos,chsize,vtitle=vtitle,vlab=vlab,time_offset=time_offset,time_scale=time_scale

; sheng. remove the default time tick. need to modify vtitle and use only vlab
if keyword_set(nouttick) then begin
    e0 = '!C'           ; enter.
    if n_elements(vtitle) eq 0 then begin
        str_element, def_opts, 'xtickformat', '(A1)', /add_replace
    endif
    if def_opts.version eq 2 and n_elements(vtitle) ne 0 then begin
        v0 = ['UT','Date']
        vs = strsplit(vtitle,e0,/extract)   ; vtitles.
        nv = n_elements(vs)
        fs = intarr(nv)                     ; flags.
        for i = 0, nv-1 do fs[i] = where(v0 eq vs[i])  ; -1: other labels.
        idx = where(fs eq -1, cnt)
        if cnt eq 0 then begin      ; no other label. xtickformat = '(A1)'
            vtitle = ''
            str_element, def_opts, 'xtickformat', '(A1)', /add_replace
        endif else begin
            vtitle = strjoin(vs[idx],e0)
            str_element, def_opts, 'xtickname', ts
            for i = 0, n_elements(ts)-1 do ts[i] = $
                strjoin((strsplit(ts[i],e0,/extract))[idx],e0)
            str_element, def_opts, 'xtickname', ts, /add_replace
        endelse
    endif
    if def_opts.version eq 3 and n_elements(vtitle) ne 0 then begin
        v0 = 'hhmm'
        vs = strsplit(vtitle,e0,/extract)   ; vtitles.
        nv = n_elements(vs)
        idx = where(vs eq v0, cnt)
        if cnt ne 0 then begin
            vs[idx:idx+1] = ''
            idx = where(vs ne '', cnt)
            if cnt eq 0 then vs = '' 
        endif else vs = ''
        
        
        if vs[0] eq '' then begin
            vtitle = ''
            str_element, def_opts, 'xtickformat', '(A1)', /add_replace
        endif else begin
            vtitle = strjoin(vs,e0)
            str_element, def_opts, 'xtickname', ts
            for i = 0, n_elements(ts)-1 do ts[i] = $
                strjoin((strsplit(ts[i],e0,/extract))[idx],e0)
            str_element, def_opts, 'xtickname', ts, /add_replace
        endelse
    endif
endif else begin    ; remove the extra line with date.
    ; sheng. remove the second line of date.
    ;if keyword_set(nodate) then begin
;    if n_elements(vtitle) eq 0 then begin
;;        str_element, def_opts, 'xtickformat', '(A1)', /add_replace
;    endif else begin
;        tpos = strpos(vtitle,'!C', /reverse_search)
;        vtitle = strmid(vtitle,0,tpos)
;    endelse
;    str_element, def_opts, 'xtickname', ts
;    tpos = strpos(ts[0],'!C', /reverse_search)
;    ts[0] = strmid(ts[0],0,tpos+2)
;    str_element, def_opts, 'xtickname', ts, /add_replace
    ;endif
endelse


; sheng.
if keyword_set(uttick) then begin
    e0 = '!C'
    if def_opts.version eq 2 and n_elements(vtitle) ne 0 then begin
        v0 = ['UT','Date']
        vs = sstrsplit(vtitle,e0,/extract)
        nv = n_elements(v0)
        fs = intarr(nv)
        for i = 0, nv-1 do fs[i] = where(vs eq v0[i])
        time_setup = time_ticks(trg)
        tickv = data_cut(uttick,time_setup.xtickv+time_offset,extrapolate=2)
        ntick = n_elements(tickv)
        get_data, uttick, tmp, val
        dt = max(val)-min(val)
        str_element, def_opts, 'xtickname', ts
        for i = 0, ntick-1 do begin
            tts = strsplit(ts[i],e0,/extract)
            tts[0] = time_string(tickv[i], tformat='hhmm')
            if i eq 0 then begin
                tts[1] = time_string(tickv[i], tformat='MTH DD')
            endif else begin
                if n_elements(tts) ge 1 then tts = [tts[0],'',tts[1:*]]
            endelse
            ts[i] = strjoin(tts,e0)
;            tmp = time_ticks(tickv[i]+[0,dt])
;            if size(tmp,/type) ne 8 then tmp = strarr(nv) else $
;            tmp = sstrsplit(tmp.xtickname[0],e0,/extract)
;            if i gt 0 then tmp[1] = ''  ; suppress date for other ticks.
;            ttick = sstrsplit(ts[i],e0,/extract)
;            for j = 0, nv-1 do ttick[j] = tmp[j]
;            ts[i] = strjoin(ttick,e0)
        endfor
        str_element, def_opts, 'xtickname', ts, /add_replace
    endif else if def_opts.version eq 3 then begin
        v0 = ['hhmm',time_string(time_offset,tformat='YYY MTH DD')]
    endif
endif



; sheng. remove this block to enable to plot versus arbitrary x.
;;return time_offset in the t_offset keyword, if requested
;if undefined(time_offset) then begin
;  dprint,'Illegal time interval.',dlevel=1  
;  return
;endif

; sheng.
if undefined(time_offset) then begin
    time_offset = 0
endif

t_offset = time_offset

def_opts.xrange = (trg-time_offset)/time_scale

if keyword_set(oplot) then def_opts.noerase = 1
if keyword_set(noerase) then def_opts.noerase = 1   ; sheng.

;for i=0,nd-1 do begin
;  polyfill,(pos[*,i])([[0,1],[2,1],[2,3],[0,3]]),color=5,/norm
;endfor

;stop

init_opts = def_opts
init_opts.xstyle = 5
;if init_opts.noerase eq 0 then erase_region,_extra=init_opts
if  init_opts.noerase eq 0 then erase
init_opts.noerase = 1
str_element,init_opts,'ystyle',5,/add
init_opts.title = ''        ; sheng, suppress title by default.


box,init_opts

def_opts.noerase = 1
str_element,tplot_vars,'options.timebar',tbnames
if keyword_set(tbnames) then begin
   tbnames = tnames(tbnames)
   ntb = n_elements(tbnames)
   for i=0,ntb-1 do begin
      t = 0
      get_data,tbnames[i],data=d
      str_element,d,'x',t
      str_element,d,'time',t
      for j=0,n_elements(t)-1 do $
         oplot,(t[j]-time_offset)/time_scale*[1,1],[0,1],linestyle=1
   endfor
endif


str_element,/add,tplot_vars,'settings.y', replicate(!y,nd)
str_element,/add,tplot_vars,'settings.clip',lonarr(6,nd)

for i=0,nd-1 do begin
   name = tplot_vars.options.varnames[i]
   def_opts.position = pos[*,i]         ;  get the correct plot position
   get_data,name,alimits=limits,ptr=pdata,data=data,index=index,dtype=dtype

   if not keyword_set(pdata) and dtype ne 3 then  dprint,verbose=verbose,'Undefined or empty variable data: ',name $
   else dprint,verbose=verbose,dlevel=1,index,name,format='(i3," ",a)'
   if keyword_set(pdata) then  nd2 = n_elements(pdata) else nd2 = 1
   if dtype eq 3 then begin
    datastr = data
    yrange = [0.,0.]
    str_element,limits,'yrange',yrange
    if ndimen(datastr) eq 0 then datastr = tnames(datastr,/all);  strsplit(datastr,/extract)
    nd2 = n_elements(datastr)
    if yrange[0] eq yrange[1] then get_ylimits,datastr,limits,trg
   endif else datastr=0

   ;allow label placing for pseudo variables
   all_labels = ''
   labflag = 0b 
   label_placement = 0 ;array to determine label positions
   labidx = 0 ;offset for indexing position array
   str_element, limits,'labflag',labflag
   if nd2 gt 1 && keyword_set(labflag) && keyword_set(datastr) then begin
     ;check for labels set on the pseudo variable, use defaults if not set
     str_element, limits,'labels',all_labels
     if ~keyword_set(all_labels) then begin
       for c=0, nd2-1 do begin
         templab = ''
         get_data, datastr[c], alimits=templim
         str_element, templim, 'labels', templab
         if keyword_set(templab) then begin
           all_labels = keyword_set(all_labels) ? [all_labels,templab]:templab
           label_placement = [label_placement,replicate(c,n_elements(templab))]
         endif
       endfor
     endif
     if n_elements(label_placement) gt 1 then begin
       label_placement = label_placement[1:n_elements(label_placement)-1]
     endif
   endif
   
   ;allow colors to be set on pseudo variables
   colors_set = 0b
   color_offset = 0
   str_element, limits, 'colors', colors_set
   
   for d=0,nd2-1 do begin
     newlim = def_opts
     newlim.ytitle = keyword_set(lazy_ytitle) ? strjoin(strsplit(name,'_',/extract),'!c')  : name
     if keyword_set(datastr) then begin
        name = datastr[d]
        get_data,name,index=index,data=data,alimits=limits2,dtype=dtype ;,ptr=pdata
        if not keyword_set(data)  then  dprint,verbose=verbose,'Unknown variable: ',name $
        else dprint,verbose=verbose,dlevel=1,index,name,format='(i3,"   ",a)'
     endif else limits2 = 0
     if size(/type,data) eq 8 then begin
        tshift = 0.d
        str_element,data,'tshift',value = tshift
        data.x = (data.x - (time_offset-tshift))/time_scale
     endif  else data={x:dindgen(2),y:findgen(2)}
     extract_tags,newlim,data,      except = ['x','y','dy','v']
     extract_tags,newlim,limits2
     extract_tags,newlim,ylimits
     extract_tags,newlim,limits
     newlim.overplot = d ne 0
     if keyword_set(overplot) then newlim.overplot = 1   ;<- *** LINE ADDED **
     if i ne (nd-1) then newlim.xtitle=''
     if i ne (nd-1) then newlim.xtickname = ' '
     
     ;add labels if set on pseudo var
     if keyword_set(all_labels) then begin
       ;labels not set on pseudo var, placement determined earlier
       if keyword_set(label_placement) then begin
         label_index = where(label_placement eq d, nl)
         if nl lt 1 then label_index = -1
       ;labels explicitly set on pseudo var, add labels in order
       endif else begin
         label_index = indgen(dimen2(data.y)) + labidx
         labidx = max(label_index) + 1
       endelse
       ;add aggregated labels/indexes for current variable
       str_element, newlim, 'label_index', label_index, /add
       str_element, newlim, 'all_labels', all_labels, /add
     endif 
     
     ;set offset into color array, if plotting pseudo vars this should
     ;allow the next variable's trace to start at the proper color
     if keyword_set(colors_set) then begin
       str_element, newlim, 'color_offset', color_offset, /add
     endif
     
     ysubtitle = struct_value(newlim,'ysubtitle',def='')
     if keyword_set(ysubtitle) then newlim.ytitle += '!c'+ysubtitle
     if newlim.spec ne 0 then routine='specplot' else routine='mplot'
     str_element,newlim,'tplot_routine',value=routine
     color_table= struct_value(newlim,'color_table',default=-1) & pct=-1
     if color_table ge 0 then loadct2,color_table,previous_ct=pct
     call_procedure,routine,data=data,limits=newlim
     if color_table ne pct then loadct2,pct
     
     ;get offset into color array (for pseudo vars)
     if keyword_set(colors_set) then begin
       str_element, newlim, 'color_offset', value=color_offset
     endif
     
   endfor
   def_opts.noerase = 1
   def_opts.title  = ''
   tplot_vars.settings.y[i]=!y
   tplot_vars.settings.clip[*,i] = !p.clip

   ; sheng, add panel labels for paper figure.
   if n_elements(figlabs) ne 0 then begin
    
       str_element,tplot_vars,'options.xcharsize', xcharsz   ; sheng
       str_element, tplot_vars, 'options.charsize', charsz
       if n_elements(xcharsz) eq 0 then xcharsz = !x.charsize
       if n_elements(charsz) eq 0 then charsz = !p.charsize
       if xcharsz eq 0 then xcharsz = 1d
       if charsz eq 0 then charsz = 1d
       vchsz = double(xcharsz*charsz)
       xspace = vchsz * !d.x_ch_size / !d.x_size
       xlmarg = def_opts.xmargin[0] & r = 2d/3
       if n_elements(vtitle) ne 0 then begin                 ; sheng
           tmp = max(strlen(strsplit(vtitle,'!C',/extract)))+2
           if tmp gt xlmarg then xlmarg = tmp           ; long vtitle.
           if tmp le xlmarg*r then xlmarg = xlmarg*r  ; short vtitle.
           if keyword_set(short_vlabel) then xlmarg = tmp  ; shortest possible.
       endif
       xpos = pos[0,nd-1] - (xlmarg-2) * xspace
       if keyword_set(vlab_marg) then xpos = pos[0,nd-1] - vlab_marg * xspace    ; set vlable pos directly
    
       if size(figlabs[i],/type) eq 8 then begin
         tlab = (figlabs[i]).text
         tcolor = (figlabs[i]).color
       endif else begin
         tlab = figlabs[i]
         tcolor = !p.color
       endelse
;       plot, [0,1],[0,1], /noerase, /nodata, xstyle = 5, ystyle = 5, $
;        position = [0,0,1,1]
       px = pos[0,i]
       py = pos[3,i]
       dx = double(!d.x_ch_size)/!d.x_size
       dy = double(!d.y_ch_size)/!d.y_size
       cx = 1.0     ; size.
       cy = 0.9
       ax = 1.0    ; padding space.
       ay = 0.1
       
       tpx = px+dx*ax
       tpy = py-dy*(cy-ay)
       
       tpx = xspace
       tpy = py-dy*0.8
;       polyfill, [px,px+dx*cx,px+dx*cx,px,px], [py,py,py-dy*cy,py-dy*cy,py], /normal, color = !p.color
       xyouts, tpx, tpy, /normal, tlab, color = tcolor, alignment = 0
;       plot, [0,1],[0,1], /noerase, /nodata, xstyle = 5, ystyle = 5, $
;        position = pos[*,i]
    endif
endfor
str_element,tplot_vars,'settings.varnames',varnames,/add_replace
str_element,tplot_vars,'settings.d',!d,/add_replace
str_element,tplot_vars,'settings.p',!p,/add_replace
str_element,tplot_vars,'settings.x',!x,/add_replace
str_element,tplot_vars,'settings.trange_cur',(!x.range * time_scale) + time_offset

;option to control left-hand labels for x-axis
str_element, def_opts, 'vtitle', vtitle

if ~keyword_set(novtitle) then begin                 ; finish var_labels
  str_element,tplot_vars,'options.xcharsize', xcharsz   ; sheng
  str_element, tplot_vars, 'options.charsize', charsz
  if n_elements(xcharsz) eq 0 then xcharsz = !x.charsize
  if n_elements(charsz) eq 0 then charsz = !p.charsize
  if xcharsz eq 0 then xcharsz = 1d
  if charsz eq 0 then charsz = 1d
  vchsz = double(xcharsz*charsz)
  xspace = vchsz * !d.x_ch_size / !d.x_size
  yspace = chsize * !d.y_ch_size / !d.y_size
  xlmarg = def_opts.xmargin[0] & r = 2d/3
  ; sheng change yyyy mth dd to yyyy and shrink to the shortest margin.
  if n_elements(vtitle) ne 0 then begin
    strs = strsplit(vtitle,'!C',/extract)
    months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    foreach tstr, strs, ii do begin
        if ii ge 1 and strs[ii-1] eq 'Seconds' then continue
        foreach month, months do begin
            index = strpos(tstr,month)
            if index eq -1 then continue
            strs[ii] = strtrim(strmid(tstr,0,index),2)
        endforeach
    endforeach
    vtitle = strjoin(strs,'!C')
    tmp = max(strlen(strs))+4
    if tmp gt xlmarg then xlmarg = tmp           ; long vtitle.
    if tmp le xlmarg*r then xlmarg = xlmarg*r  ; short vtitle.
    xlmarg = tmp  ; shortest possible.
  endif
  xpos = pos[0,nd-1] - xlmarg * xspace
  if keyword_set(vlab_marg) then xpos = pos[0,nd-1] - vlab_marg * xspace    ; set vlable pos directly
  name = tplot_vars.options.varnames[nd-1]
  get_data, name, limits=lim
  str_element,lim,'xticklen', value=xticklen
  if n_elements(xticklen) eq 0 then str_element, tplot_vars.options, 'xticklen', value=xticklen
  if n_elements(xticklen) eq 0 then str_element,tplot_vars,'settings.x', value=xticklen
  if n_elements(xticklen) eq 0 then xticklen = 0
  if xticklen ge 0 then xticklen = 0
  ypos = pos[1,nd-1] - 1.5 * yspace + xticklen*(pos[3,nd-1]-pos[1,nd-1])
  if n_elements(vtitle) eq 0 then vtitle = ''
  xyouts,xpos,ypos,vtitle,/norm,charsize=vchsz
endif

;time_stamp,charsize = chsize*.5

if (!d.flags and 256) ne 0  then begin    ; windowing devices
  str_element,tplot_vars,'settings.window',!d.window,/add_replace
  if def_opts.window ge 0 then wset,current_window
endif
!x = plt.x
!y = plt.y
!z = plt.z
!p = plt.p


str_element,tplot_vars,'settings.time_scale',time_scale,/add_replace
str_element,tplot_vars,'settings.time_offset',time_offset,/add_replace
new_tvars = tplot_vars
return
end
