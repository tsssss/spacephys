
utr = time_double(['2010-01-01','2017-01-01'])
ut0 = 100d*86400*365

; load sunspot number.
load = 0
if tnames('ssn') eq '' then begin
    ssnfn = srootdir()+'/SN_d_tot_V2.0.txt'
    nline = file_lines(ssnfn)
    lines = strarr(nline)
    openr, lun, ssnfn, /get_lun
    readf, lun, lines
    free_lun, lun
    
    
    uts = dblarr(nline)
    ssns = dblarr(nline)
    dssns = dblarr(nline)
    yrs = dblarr(nline)
    for i = 0, nline-1 do begin
        tline = lines[i]
        tmp = strsplit(tline,' ', /extract)
        if strlen(tmp[1]) lt 2 then tmp[1] = '0'+tmp[1]
        if strlen(tmp[2]) lt 2 then tmp[2] = '0'+tmp[2]
        uts[i] = time_double(tmp[0]+tmp[1]+tmp[2],tformat='YYYYMMDD')
        yrs[i] = float(tmp[3])
        ssns[i] = float(tmp[4])
        dssns[i] = float(tmp[5])
    endfor
    
    idx = where(ssns eq -1, cnt)
    if cnt ne 0 then begin
        ssns[idx] = !values.f_nan
        dssns[idx] = !values.f_nan
    endif
    
    idx = where(uts ge utr[0] and uts le utr[1])
    store_data, 'ssn', uts, ssns, dssns, limits = $
        {ytitle:'Sun Spot Num'}
endif



; **** calc MAT, filter, etc.


idx = where(uts ge utr[0] and uts le utr[1])
x0 = yrs[idx]
f0 = ssns[idx]


zr = [-1,1]*15
nlevel = 20
levels = smkarthm(zr[0],zr[1],nlevel,'n')
colors = smkarthm(10,250,nlevel,'n')
yr = [4,400]
xr = minmax(x0)
order = 2

dr0 = 86400d     ; 1 day.
sclinfo = [1,500,40]    ; 5 days to more than 11 yr.
rscales = smkgmtrc(sclinfo[0],sclinfo[1],sclinfo[2],'n')
rscales = suniq(floor(rscales))
tscales = rscales*dr0

f0 = smooth(f0,sclinfo[0],/nan)
f0mat = swvmat(f0, order, scale = rscales)

ofn = shomedir()+'/test_mat_ssn.pdf'
ofn = 0
sgopen, ofn, xsize = 12, ysize = 8, /inch
device, decomposed = 0
loadct, 39

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size


erase, 255

poss = sgcalcpos(2, ypad = 8)

plot, x0, f0, /noerase, position = poss[*,0], $
    xstyle = 1, xrange = xr, xtickformat = '(A1)', $
    ystyle = 1, ytitle = 'Sun Spot Num', color = 0

contour, f0mat, x0, rscales, /fill, position = poss[*,1], /noerase, $
    /ylog, ystyle = 1, yrange = yr, ytitle = 'Period!C(Day)', $
    xstyle = 1, xrange = xr, xtitle = 'Time (Year)', $
    levels = levels, nlevel = nlevel, c_colors = colors, background = 0

lines = [365,60,27,27*0.5]
foreach tline, lines do begin
    plots, xr, tline+[0,0], color = 255
    tx = xr[1] & ty = tline
    tmp = convert_coord(tx,ty, /data, /to_normal)
    xyouts, tmp[0]+xchsz*1, tmp[1], /normal, sgnum2str(tline)+' days', color = 0
endforeach

tpos = poss[*,1] & tpos[1] = tpos[3]-ychsz*1 & tpos = tpos+[0,1,0,1]*ychsz*2
sgcolorbar, colors, /horizontal, zrange = zr, position = tpos, $
    ztitle = 'Sun Spot Num', zcharsize = 1
    
sgclose


end