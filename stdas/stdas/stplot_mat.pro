;+
; Type: procedure.
; Purpose: Wrap swvmat for quick spectrogram within tplot.
; Parameters: vname, in, type = string, required. Signal varname in tplot.
; Keywords: scale, in/out, type = [m], optional. Scale in time.
;   nscale, in/out, integer, optional. Set # of scales.
;   order, in, integer, optional. Order of mat, see swvmat.
;   newname, in, string, optional. Output varname in tplot.
;   overwrite, in, boolean, optional. Set to overwrite vname.
;   ytitle, in, stirng, optional. Title for var.
; Notes: The variable should be 1-d data. This code forces label to be
;   specified by label keyword if vname contains no label.
; Dependence: tplot,slib.
; Author: Sheng Tian.
; History: 2013-11-01, Sheng Tian, create.
;-
pro stplot_mat, vname, scale = tscales, nscale = nscale, $
    order = order, newname = newname, overwrite = overwrite, $
    ytitle = ytitle, zrange = zrange

    get_data, vname, t0, f0, limits = lim
    dr = sdatarate(t0) & dr1 = 1d/dr
    nrec = n_elements(t0)
    
    if n_elements(tscales) eq 0 then begin
        nscale = keyword_set(nscale)? nscale: 60
        rscales = round(smkgmtrc(4,0.5*nrec,nscale,'n'))
        rscales = double(rscales[uniq(rscales)])
        tscales = rscales*dr
    endif else rscales = tscales*dr1
    nscale = n_elements(tscales)
    yrange = minmax(tscales)

    if keyword_set(overwrite) then newname = vname
    if ~keyword_set(newname) then newname = vname+'_mat'
    if ~keyword_set(ytitle) then ytitle = 'MAT (s)'
    if ~keyword_set(zrange) then zrange = [min(f0),max(f0)]/nscale*5
    ztitle = stagexist(lim,'ytitle')? lim.ytitle: ''

    f0mat = swvmat(f0, order, scale = rscales)
    tscales = rscales*dr
    
    store_data, newname, t0, f0mat, tscales, $
        limits = {spec:1, no_interp:1, ytitle:ytitle, zrange:zrange}
    ylim, newname, yrange[0], yrange[1], 1

end