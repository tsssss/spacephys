
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
perp = '!9'+string(94b)+'!X'
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




; **** plot 1: EB1 and MB1 in FAC.

    ofn = 0
    ofn = shomedir()+'/es_wave/fig_fields_fac.pdf'
    sgopen, ofn, xsize = 5, ysize = 4, /inch
    
    vars = pre0+['eb1','mb1']+'_fac'
    nvar = n_elements(vars)
    
    poss = sgcalcpos(nvar, lmargin=15, rmargin=5)
    xchsz = double(!d.x_ch_size)/!d.x_size
    ychsz = double(!d.y_ch_size)/!d.y_size
    
    options, vars, 'labels', ['para','west','north']
    tvar = pre0+'eb1_fac'
    options, tvar, 'yrange', [-1,1]*110
    options, tvar, 'yticks', 2
    options, tvar, 'ytickv', [-1,0,1]*100
    options, tvar, 'yminor', 5
    
    tvar = pre0+'mb1_fac'
    options, tvar, 'yrange', [-1,1]*0.13
    options, tvar, 'yticks', 2
    options, tvar, 'ytickv', [-1,0,1]*0.1
    options, tvar, 'yminor', 5

    
    get_data, vars[0], uts
    tplot, vars, trange = minmax(uts), position = poss
    
    figlabs = ['a. dE FAC','b. dB FAC']
    for i = 0, nvar-1 do xyouts, $
        poss[0,i]-xchsz*12, poss[3,i]-ychsz*1, /normal, figlabs[i]

    sgclose



end
