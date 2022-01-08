
faclabs = ['||','!9^!X,West','!9^!X,North']

ifn = '/Users/Sheng/Google Drive/works/data/2013_0607/double_layer_flag.csv'
nline = file_lines(ifn)
lines = strarr(nline)

openr, lun, ifn, /get_lun
readf, lun, lines
free_lun, lun

nrec = nline*0.5

uts = dblarr(nrec)
flags = dblarr(nrec)
for i = 0, nrec-1 do begin
    tmp = strsplit(lines[i],',',/extract)
    uts[i] = time_double(strmid(tmp[0],24,20),tformat='YYYY_MMDD_hhmm_ss.ff')
    flags[i] = fix(tmp[1])
endfor
idx = where(flags eq 1, cnt)
if cnt ne 0 then flags[idx] = 0
store_data, 'rbspa_dblayer_flag', uts, flags, limits = $
    {yrange:[-1,3], ystyle:1, ytitle:'Double Layer!C2=exist', ytickv:[0,2], yticks:2}


for i = 0, nrec-1 do begin
    tmp = strsplit(lines[i+nrec],',',/extract)
    uts[i] = time_double(strmid(tmp[0],24,20),tformat='YYYY_MMDD_hhmm_ss.ff')
    flags[i] = fix(tmp[1])
endfor
idx = where(flags eq 1, cnt)
if cnt ne 0 then flags[idx] = 0
store_data, 'rbspb_dblayer_flag', uts, flags, limits = $
    {yrange:[-1,3], ystyle:1, ytitle:'Double Layer!C2=exist', ytickv:[0,2], yticks:2}


options, 'rbsp?_de_fac_mat1', 'ytitle', 'dE FAC!C1/16-4 sec!C(mV/m)'
options, 'rbsp?_de_fac_mat1', 'labels', 'dE!D'+faclabs

tvar = ['vsc','dblayer_flag','de_fac_mat1']
nvar = n_elements(tvar)

options, 'rbspa_vsc', 'yrange', [-30,0]

ofn = shomedir()+'/fig_dblayer_efield.pdf'
;ofn = 0
sgopen, ofn, xsize = 6, ysize = 8, /inch

device, decomposed = 0
loadct2, 43

xchsz = double(!d.x_ch_size)/!d.x_size
ychsz = double(!d.y_ch_size)/!d.y_size


poss = sgcalcpos(nvar, region=[0,0.5,1,1])
vars = 'rbspa_'+tvar
tplot, vars, position = poss, /noerase, /novtitle
xyouts, 0.5*(poss[0,0]+poss[2,0]), poss[3,0]+ychsz*0.2, /normal, alignment = 0.5, 'RBSP-A'



poss = sgcalcpos(nvar, region=[0,0,1,0.5])
vars = 'rbspb_'+tvar
tplot, vars, position = poss, /noerase, /novtitle
xyouts, 0.5*(poss[0,0]+poss[2,0]), poss[3,0]+ychsz*0.2, /normal, alignment = 0.5, 'RBSP-B'

sgclose

end
