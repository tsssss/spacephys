;+
; Type: procedure.
; Purpose: Calculate the moving average transform and the wave bands.
; Parameters: vname, in, string, req. Signal varname in tplot, must be 1d.
; Keywords:
;   order, in, integer, opt. Order of mat, see swvmat. Default order 2.
;   scale, in/out, dblarr[m], opt. Scale in time.
;   scaleinfo, in/out, double/dblarr[3], opt. Set nscale or [min,max,nscale].
; Notes: The variable should be 1-d data.
; Dependence: tplot,slib.
; History:
;   2014-11-05, Sheng Tian, create.
;-
pro stplot_mat_filter, vname, order = order, $
    scale = tscs, scaleinfo = tscinfo, filter = tfts, zrange = zrange

    ; get original data.
    get_data, vname, t0, f0, limits = lim
    dr = sdatarate(t0) & dr1 = 1d/dr
    nrec = n_elements(t0)
    tr = t0[[0,nrec-1]]
    
    ; plot settings.
    device, decomposed = 0
    loadct2, 43
    
    if n_elements(matname) eq 0 then matname = vname+'_mat'
    
    ; flag for interactively determining time scale info, zrange, filters.
    nudge = 'n'
    
    ; check mat scale, want tscs, mintsc, maxtsc, nsc.
    if n_elements(tscs) ne 0 then begin             ; have time scale.
        mintsc = min(tscs, max = maxtsc)
        nsc = n_elements(tscs)
    endif else begin
        mintsc = 4*dr & maxtsc = 0.5*nrec*dr & nsc = 50
        case n_elements(tscinfo) of
            3: begin                                ; full time scale info.
                mintsc = tscinfo[0] & maxtsc = tscinfo[1]
                nsc = tscinfo[2] & end
            2: begin                                ; min,max scale in time.
                mintsc = tscinfo[0] & maxtsc = tscinfo[1] & end
            1: begin                                ; nscale.
                nsc = tscinfo[2] & end
            0: nudge = 'y'
        endcase
        tscs = smkgmtrc(mintsc,maxtsc,nsc,'n')
    endelse
    
    ; check mat zrange.
    if n_elements(zrange) eq 0 then begin
        nudge = 'y'
        zrange =  [min(f0),max(f0)]/nsc
    endif
    
    ; do mat.
    stplot_mat, vname, newname = matname, scale = tscs, zrange = zrange
    mintsc = min(tscs, max = maxtsc) & nsc = n_elements(tscs)
    
    ; determine tscs, mintsc, maxtsc, nsc, zrange.
    if nudge eq 'y' then begin
        vars = [vname,matname]
        tplot, vars
        while nudge eq 'y' do begin
            printf, -1, 'Current scale info: '
            printf, -1, [mintsc,maxtsc,nsc]
            tmp = ''
            read, tmp, prompt = 'Time scale info? (min, max, n): '
            if tmp ne '' then begin
                tscinfo = double(strsplit(tmp,', ',/extract))
                case n_elements(tscinfo) of
                    3: begin                            ; full time scale info.
                        mintsc = tscinfo[0] & maxtsc = tscinfo[1]
                        nsc = tscinfo[2] & end
                    2: begin                            ; min,max scale in time.
                        mintsc = tscinfo[0] & maxtsc = tscinfo[1] & end
                    1: begin                            ; nscale.
                        nsc = tscinfo[2] & end
                endcase
                tscs = smkgmtrc(mintsc,maxtsc,nsc,'n')
                stplot_mat, vname, newname = matname, scale = tscs
                mintsc = min(tscs, max = maxtsc) & nsc = n_elements(tscs)
            endif
            ; read zrange.
            printf, -1, 'Current zrange: '
            printf, -1, zrange
            tmp = ''
            read, tmp, prompt = 'Zrange? (min,max): '
            if tmp ne '' then zrange = double(strsplit(tmp,', ',/extract))
            options, matname, 'zrange', zrange
            tplot, vars
            read, nudge, prompt = 'Tune up scale? (y/n): '
            if nudge eq '' then nudge = 'y'
        endwhile
    endif
    
    ; check filters. want fts, nft, ftids.
    if n_elements(tfts) le 1 then begin
        nudge = 'y' & tfts = [mintsc,maxtsc]
    endif
    tfts = stofilter(tfts) & nft = n_elements(tfts)/2
    ftids = string(indgen(nft)+1,format='(I0)')
    
    ; do filter.
    for i = 0, nft-1 do $
        stplot_filter, matname, filter = tfts[i,*], ifilter = ftids[i]
    
    if nudge eq 'y' then begin
        vars = [matname,matname+'f'+ftids]
        tplot, vars
        while nudge eq 'y' do begin
            store_data, matname+'f'+ftids, /delete
            ctime, tmp, tfts, /exact
            if n_elements(tfts) ge 2 then begin
                store_data, matname+'f'+ftids, /delete
                tfts = stofilter(tfts) & nft = n_elements(tfts)/2
                ftids = string(indgen(nft)+1,format='(I0)')
                for i = 0, nft-1 do $
                    stplot_filter, matname, filter = tfts[i,*], ifilter = ftids[i]
                vars = [matname,matname+'f'+ftids]
                tplot, vars
            endif
            printf, -1, 'Current filter: '
            printf, -1, [tfts[*,0],tfts[nft*2-1]]
            read, nudge, prompt = 'Tune up filter? (y/n): '
            if nudge eq '' then nudge = 'y'
        endwhile
    endif
    
    ; more on filters.
    store_data, matname+'_filter', tr, [tfts[*,0],tfts[nft*2-1]]##[1,1], $
        limits = {colors:intarr(nft+1)-1}
    store_data, matname+'_comb', data = matname+['','_filter'], $
        limits = {yrange:[mintsc,maxtsc],ylog:1}
    ; leftover wave.
    stplot_filter, matname, filter = [tfts[0],tfts[nft*2-1]], ifilter = 'del'
    get_data, matname+'fdel', t0, f1
    get_data, vname, t0, f0
    store_data, matname+'fdel', t0, f0-f1
    
    ; plot all data.
    vars = [matname+'_comb', matname+'f'+ftids, matname+'fdel', vname]
    tplot, vars
    printf, -1, 'Current scale info: '
    printf, -1, [mintsc,maxtsc,nsc]
    printf, -1, 'Current zrange: '
    printf, -1, zrange
    printf, -1, 'Current filter: '
    printf, -1, [tfts[*,0],tfts[nft*2-1]]
    
end


fn = shomedir()+'/Dropbox/code/slib/topics/qshi/miro/xx++.itx'
nrec = file_lines(fn)
header = strarr(3)
openr, lun, fn, /get_lun
readf, lun, header
tmp = dblarr(2,nrec-6)
readf, lun, tmp
free_lun, lun
x0 = reform(tmp[0,*])
y0 = reform(tmp[1,*])

nrec = n_elements(x0)
dr = (max(x0)-min(x0))/(nrec-1)
x1 = smkarthm(min(x0),max(x0),nrec,'n')
dr = sdatarate(x1)
y1 = sinterpol(y0,x0,x1)

store_data, 'rho', x1, y1, limits = {yrange:[min(y1),max(y1)]}
stplot_mat_filter, 'rho', scaleinfo = [0.001,0.05,50], zrange = 0.01*[-1,1]
end
