
pro miro_qshi, fn

    nrec = file_lines(fn)
    openr, lun, fn, /get_lun
    tmp = dblarr(2,nrec) & readf, lun, tmp
    free_lun, lun
    y0 = reform(tmp[0,*])
    x0 = reform(tmp[1,*])

    nrec = n_elements(x0)
    dr = (max(x0)-min(x0))/(nrec-1)
    x1 = smkarthm(min(x0),max(x0),nrec,'n')
    dr = sdatarate(x1)
    y1 = sinterpol(y0,x0,x1)
    
    minsc = 0.5    ; 6 waves over 10.
    maxsc = 15
    nsc = 50
    dat_sc = smkgmtrc(minsc, maxsc, nsc, 'n')
    rec_sc = dat_sc/dr
    x1_mat = swvmat(y1, scale = rec_sc)
    dat_sc = rec_sc*dr; & dat_sc = findgen(n_elements(xx)) # dat_sc
    
    ; mat.
    store_data, 'rho', x1, y1, limits = {spec:0, yrange:[min(y1),max(y1)]}
    store_data, 'rho_mat', x1, x1_mat, dat_sc, limits = $
        {spec:1, yrange:[min(dat_sc),max(dat_sc)], ystyle:1, ylog:1, zrange:[-3,3]}
    
    ; filter.
    fltr = dblarr(nrec)+2
    miro = dblarr(nrec)
    sdho = dblarr(nrec)
    for i = 0, nrec-1 do begin
        idx = (where(dat_sc gt fltr[i]))[0]
        sdho[i] = total(x1_mat[i,0:idx])
        miro[i] = total(x1_mat[i,idx:*])
    endfor
    store_data, 'miro', x1, miro, limits = {yrange:[-40,40],ysytle:1}
    store_data, 'sdho', x1, sdho, limits = {yrange:[-40,40],ysytle:1}
    bgrho = y1-miro-sdho
    store_data, 'background', x1, bgrho
    store_data, 'combine', x1, [[miro+bgrho],[sdho]], limits = $
        {colors:[6,4],labels:['miro+bg','sdho'],yrange:[-40,40],ystyle:1}
        
    device, decomposed = 0
    loadct2, 43
    !p.font = 1
    tplot_options, 'ygap', 0.25
    tplot_options, 'ynozero', 1
    tplot_options, 'version', 2
    tplot_options, 'num_lab_min', 8
    tplot_options, 'labflag', 1
    tplot_options, 'zcharsize', 0.9
    time_stamp, /off
    tplot, ['rho','rho_mat','miro','sdho','background','combine']
    pstplot, filename = '~/miro_qshi.eps'
    
    ; output, rho, miro, sdho.
    datfn = '~/miro_tian.dat'
    openw, lun, datfn, /get_lun
    printf, lun, 'x y miro sdho bg'
    for i = 0, nrec-1 do begin
        printf, lun, x1[i], y1[i], miro[i], sdho[i], bgrho[i]
    endfor
    free_lun, lun
end

fn = srootdir()+'QS.dat'
miro_qshi, fn
end