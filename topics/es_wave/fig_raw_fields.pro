
; settings.
utr0 = time_double('2013-06-07/04:56'+[':02.00',':03.00'])  ; total data.
utr1 = time_double('2013-06-07/04:56'+[':02.20',':02.80'])  ; zoom in for wave packet.


tprobe='a'
pre0 = 'rbsp'+tprobe+'_'
dr0 = 1d/4096
tpad = 200*dr0  ; sec.
p0 = 0.02       ; sec, periods of the main wave.
utr = utr0+[-1,1]*tpad
timespan, utr[0], utr[1]-utr[0], /second

re = 6378d & re1 = 1d/re
times = '!9'+string(180b)+'!X'
re = 6378d & re1 = 1d/re
rgb = sgcolor(['red','green','blue'])
c6s = sgcolor(['magenta','blue','cyan','green','yellow','red'])
uvw = ['U','V','W']
xyz = ['X','Y','Z']


tplot_options, 'constant', 0
tplot_options, 'labflag', -1
tplot_options, 'num_lab_min', 5
tplot_options, 'version', 1
tplot_options, 'ystyle', 1

; rbspx_[vb1,eb1,mb1]. original E/B field in UVW.
; rbspx_[uvw]_gsm. UVW direction in GSM.
if tnames(pre0+'vb1') eq 0 then _2013_0607_0456_load_burst_data


vars = pre0+'vb1'
vars = pre0+['eb1','mb1']

vars = pre0+['eb1','mb1']+'_uvw'


tvar = pre0+'vb1'
get_data, tvar, uts, vb1
vsc = [[vb1[*,0]+vb1[*,1]],[vb1[*,2]+vb1[*,3]],[vb1[*,4]+vb1[*,5]]]*0.5
store_data, pre0+'vsc', uts, vsc, limits = $
    {ytitle:'(V)', labels:'Vsc '+uvw, colors:rgb}


fit_uv = linfit(vsc[*,0],vsc[*,1])  ; -0.05 offset, 1.03 slope.
fit_uw = linfit(vsc[*,0],vsc[*,2])  ; -3.55 offset, 0.40 slope.

for i = 0, 2 do begin
    vsc[*,i] = vsc[*,i]-vsc[0,i]
    if i eq 2 then vsc[*,i] = vsc[*,i]*2.5
endfor
store_data, pre0+'vsc1', uts, vsc, limits = $
    {ytitle:'(V)', labels:'Vsc '+uvw+['','','!C  '+times+'2.5'], colors:rgb}





 ; **** plot 2: raw vsc, eb1, vb1, b0.
    b0w = 0.2   ; width to smooth over.
    get_data, pre0+'b0_gsm', uts, dat
    b0 = snorm(dat)
    b1 = smooth(b0,b0w/sdatarate(uts),/edge_mirror)
    store_data, pre0+'b0', uts, [[b0],[b1]], limits = $
        {ytitle:'(nT)', labels:['|B|','smth.'], ynozero:1, yticks:3, colors:sgcolor(['black','red'])}
    vars = pre0+['vsc1','eb1_uvw','mb1_uvw','b0']
    nvar = n_elements(vars)
    
    ofn = 0
    ofn = shomedir()+'/es_wave/fig_raw_vsc_eb_b0.pdf'
    sgopen, ofn, xsize = 5, ysize = 6, /inch
    
    poss = sgcalcpos(nvar, lmargin=15, rmargin=5)
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    tvar = pre0+'vsc1'
    options, tvar, 'yrange', [-7,1]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 4
    
    tvar = pre0+'eb1_uvw'
    options, tvar, 'yrange', [-1,1]*110
    options, tvar, 'yticks', 2
    options, tvar, 'ytickv', [-1,0,1]*100
    options, tvar, 'yminor', 5

    tvar = pre0+'mb1_uvw'
    options, tvar, 'yrange', [-1,1]*0.13
    options, tvar, 'yticks', 2
    options, tvar, 'ytickv', [-1,0,1]*0.1
    options, tvar, 'yminor', 5
    
    tvar = pre0+'b0'
    options, tvar, 'yrange', [240.2,241.8]
    options, tvar, 'yticks', 2
    options, tvar, 'yminor', 4
    
    get_data, vars[0], uts
    tplot, vars, trange = minmax(uts), position = poss
    
    figlabs = ['a. V!DSC','b. dE UVW','c. dB UVW','d. |B|']
    for i = 0, nvar-1 do xyouts, $
        poss[0,i]-xchsz*12, poss[3,i]-ychsz*1, /normal, figlabs[i]
        
    sgclose
    



; **** plot 1: raw VB1.
    tvar = pre0+'vb1'
    get_data, tvar, uts, vb1
    for i = 0, 5 do begin
        id = string(i+1,format='(I0)')
        store_data, tvar+id, uts, vb1[*,i], limits = $
            {ytitle:'(V)', labels:'V'+id}
    endfor
    
    ofn = 0
    ofn = shomedir()+'/es_wave/fig_raw_vb1.pdf'
    sgopen, ofn, xsize = 5, ysize = 8, /inch
    
    poss = sgcalcpos(6, lmargin=15, rmargin=5)
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    
    vars = tnames(tvar+'?')
    nvar = n_elements(vars)
    options, vars, 'yticks', 2
    options, vars, 'ystyle', 1
    options, vars, 'labels', ''
    options, vars, 'yminor', 5
    options, vars[0:1], 'yrange', [-24,-9]
    options, vars[2:3], 'yrange', [-24,-9]
    options, vars[4:5], 'yrange', [-12,-7]
    tplot, vars, trange = minmax(uts), position = poss
    
    figlabs = ['a','b','c','d','e','f']+'.'
    for i = 0, nvar-1 do xyouts, $
        poss[0,i]-xchsz*10, poss[3,i]-ychsz*1, /normal, figlabs[i]+' V'+string(i+1,format='(I0)')
    
    sgclose



end
