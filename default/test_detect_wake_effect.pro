;+
; Plot waveform and morlet wavelet for a whole day, to show overall features in time-frequency domain.
;-


;---Load data.
test = 0
secofday = 86400d

date = '1998-09-25'
;date = '1998-09-26'
;date = '1998-09-14'
mission = 'polar'
probe = ''
var = 'po_e_uv'

date = '2013-05-01'
;date = '2013-06-07'
probe = 'b'
;    date = '2012-11-14'
;    probe = 'a'
mission = 'rbsp'
var = 'rbsp'+probe+'_e_uv'

case mission of
    'polar': spin_rate = polar_info('spin_rate')
    'rbsp': spin_rate = rbsp_info('spin_rate')
    'themis': spin_rate = themis_info('spin_rate')
endcase

duration = 10*spin_rate
full_time_range = time_double(date)+[0,secofday]
nsection = round((full_time_range[1]-full_time_range[0])/duration)
sections = smkarthm(full_time_range[0], duration, nsection, 'x0')


spin_freq = 1d/spin_rate
base = 2d^(1d/4)
freq_range = [0.5,8]/spin_rate
selected_freqs = smkgmtrc(min(freq_range),max(freq_range),base,'dx')
nselected_freq = n_elements(selected_freqs)
selected_periods = 1d/selected_freqs
wavelet_info = wavelet_info()
t2s = wavelet_info[7]
selected_scales = selected_periods*t2s
if tnames(var) ne '' then begin
    get_data, var, times
    index = where_pro(times,full_time_range, count=count)
    if count le 1 then store_data, var, /delete
endif
if tnames(var) eq '' then begin
    case mission of
        'polar': polar_read_euv, full_time_range
        'rbsp': rbsp_read_euv, full_time_range, probe=probe
        'themis': themis_read_euv, full_time_range, probe=probe
    endcase
endif

;---Do wavelet analysis at given frequencies.
u_var = var+'_section'
u_index = 0
cwt_var = u_var+'_cwt'
if tnames(u_var) ne '' then begin
    get_data, u_var, uts
    index = where_pro(uts, full_time_range, count=count)
    if count lt nsection then store_data, u_var, /delete
endif
if tnames(u_var) eq '' then begin
    get_data, var, times, data
    stplot_index, var, u_index, newnames=u_var
    get_data, u_var, times, data
    yrange = [-1,1]*max(abs(sg_autolim(data)))
    add_setting, u_var, /smart, {$
        yrange: yrange, $
        short_name: 'UVW E!Du!N', $
        display_type: 'scalar'}
    calc_psd, u_var, scales=selected_scales
endif
if tnames(cwt_var) ne '' then begin
    get_data, cwt_var, 0, cwt
    if n_elements(cwt.s_j) ne nselected_freq then calc_psd, u_var, scales=selected_scales
endif


;---Plot settings.
sgopen, 0, xsize=4, ysize=3, magnify=magnify
poss = sgcalcpos(1, tmargin=2, bmargin=4, rmargin=2, lmargin=8, xchsz=xchsz, ychsz=ychsz)
sgclose, /wdelete
tpos = poss[*,0]

white = 255
red = 6
const_periods = spin_rate/[1.,2,3,4]
xticklen = -0.06
yticklen = -0.01
unit = 'mV/m'
xtitle = 'Period (sec)'
xticklen = xticklen
xlog = 1
ytitle = '|'+unit+'|/Hz!U1/2!N'
yticklen = yticklen
ylog = 1
yrange = [-3,3]
ytickv = make_bins(yrange,1)
yrange = 10.^minmax(ytickv)
yticks = n_elements(ytickv)-1
ytickn = strarr(yticks+1)
for ii=0, yticks do begin
    case ytickv[ii] of
        0: ytickn[ii] = '1'
        else: ytickn[ii] = '10!U'+string(ytickv[ii],format='(I0)')
    endcase
endfor
if yticks ge 6 then begin
    target = 3
    for ii=0, yticks do begin
        if ii mod target eq 0 then continue
        ytickn[ii] = ' '
    endfor
endif
ytickv = 10.^ytickv
yminor = 10

xrange = [0.5,20]
xtickv = [1.,10]

label_size = 0.7
psym = 3    ; dot.
symsize = 0.1




;---Calculate the amplitude.
get_data, u_var, times, edata
get_data, cwt_var, 0, cwt
ww = cwt.w_nj
s2t = cwt.s2t
cdelta = cwt.cdelta
dt = cwt.dt
amps = dblarr(nsection,nselected_freq)
foreach time, sections, ii do begin
    section_time_range = time+[0,duration]
    index = where_pro(times, section_time_range, count=N)
    uts = times[index]
    ees = edata[index]

    wps = abs(ww[index,*])^2
    gws = total(wps,1)/N
    psd = gws*s2t*dt/cdelta
    amps[ii,*] = sqrt(psd)
endforeach
amp_var = var+'_amp'
store_data, amp_var, sections, amps, selected_periods
add_setting, amp_var, /smart, {$
    ytitle: 'Period (sec)', $
    yrange: xrange, $
    zrange: [0.5,50], $
    zlog: ylog, $
    display_type: 'spec', $
    unit: '['+unit+'] /Hz!U1/2!N', $
    constant: const_periods, $
    const_color: white, $
    short_name: '|E|'}
options, amp_var, 'no_interp', 0
options, amp_var, 'x_no_interp', 1
options, amp_var, 'y_no_interp', 0

;---Check peak at f, i.e., at the spin rate.
freq_1f = min(selected_freqs-spin_freq, /absolute, index_1f)
has_1f_peak = amps[*,index_1f] ge amps[*,index_1f-2] and amps[*,index_1f] ge amps[*,index_1f+2]
flag_1f_var = var+'_has_1f_peak'
store_data, flag_1f_var, sections, has_1f_peak
add_setting, flag_1f_var, /smart, {$
    ytitle: '', $
    yrange: [-0.2,1.2], $
    yticks: 1, $
    ytickv: [0,1], $
    panel_size: 0.5, $
    display_type: 'scalar', $
    short_name: 'peak@f'}

;---Check peak at 2f.
freq_2f = min(selected_freqs-2*spin_freq, /absolute, index_2f)
has_2f_peak = amps[*,index_2f] ge amps[*,index_2f-1] and amps[*,index_2f] ge amps[*,index_2f+1]
flag_2f_var = var+'_has_2f_peak'
store_data, flag_2f_var, sections, has_2f_peak
add_setting, flag_2f_var, /smart, {$
    ytitle: '', $
    yrange: [-0.2,1.2], $
    yticks: 1, $
    ytickv: [0,1], $
    panel_size: 0.5, $
    display_type: 'scalar', $
    short_name: 'peak@2f'}

;---Check peak at 3f
freq_3f = 3*spin_freq
amp_3f = dblarr(nsection)
for ii=0, nsection-1 do begin
    amp_3f[ii] = interpol(amps[ii,*],selected_freqs, freq_3f, /quadratic)
endfor
tmp = min(selected_freqs-2.5*spin_freq, /absolute, index_3fm)
amp_3fm = amps[*,index_3fm]
tmp = min(selected_freqs-3.5*spin_freq, /absolute, index_3fp)
amp_3fp = amps[*,index_3fp]
has_3f_peak = amp_3f ge amp_3fm and amp_3f ge amp_3fp
flag_3f_var = var+'_has_3f_peak'
store_data, flag_3f_var, sections, has_3f_peak
add_setting, flag_3f_var, /smart, {$
    ytitle: '', $
    yrange: [-0.2,1.2], $
    yticks: 1, $
    ytickv: [0,1], $
    panel_size: 0.5, $
    display_type: 'scalar', $
    short_name: 'peak@3f'}

;---Check the height of the peak at 1 and 3 f.
file = join_path([srootdir(),'plot','fig_eu_wavelet_'+mission+probe+'_'+time_string(full_time_range[0],tformat='YYYY_MMDD')+'.pdf'])
if keyword_set(test) then file = test
sgopen, file, xsize=10, ysize=4
device, decomposed=0
loadct2, 43

vars = [u_var,amp_var]
nvar = n_elements(vars)
poss = sgcalcpos(nvar, ypans=[2,1.5], rmargin=8, tmargin=1)

tvar = u_var
ypans = [1,2,1]
lin_yrange = [-1,1]*10
lin_ytickv = [-2,-1,0,1,2]*5
lin_yticks = n_elements(lin_ytickv)-1
lin_ytickn = string(lin_ytickv,format='(I0)')
lin_yminor = 5

log_yrange = [-1,1]*1000
log_ytickv = 10.^[1,2,3]
log_yticks = n_elements(log_ytickv)-1
log_ytickn = string(log_ytickv,format='(I0)')
log_yminor = 10

tpos = poss[*,0]
tx = xchsz*0.5
ty = tpos[3]-ychsz*0.5
xyouts, tx,ty,/normal, 'a. E!Du'
tposs = sgcalcpos(3, ypans=ypans, ypad=0, position=tpos)
get_data, tvar, xx, yy, limits=lim

; negative part.
tpos = tposs[*,2]
yrange = -[min(log_yrange),min(lin_yrange)]
ylog = 1
ytickn = '-'+log_ytickn & ytickn[0] = ' '
store_data, tvar+'_neg', xx, -yy, limits=lim
store_data, tvar+'_neg', limits={$
    ytitle: '', $
    labels: '', $
    xticklen: xticklen, $
    ytickv: log_ytickv, $
    yticks: log_yticks, $
    yminor: log_yminor, $
    ytickname: ytickn, $
    yrange: yrange, $
    yticklen: yticklen, $
    ylog: ylog}
tplot, tvar+'_neg', trange=full_time_range, position=tpos, /noerase, /nouttick
store_data, tvar+'_neg', /delete

; posivie part.
tpos = tposs[*,0]
yrange = [max(lin_yrange),max(log_yrange)]
ylog = 1
ytickn = log_ytickn & ytickn[0] = ' '
store_data, tvar, limits={$
    ytitle: '', $
    labels: '', $
    xticklen: xticklen, $
    ytickv: log_ytickv, $
    yticks: log_yticks, $
    yminor: log_yminor, $
    ytickname: ytickn, $
    yrange: yrange, $
    yticklen: yticklen, $
    ylog: ylog}
tplot, tvar, trange=full_time_range, position=tpos, /noerase, /nouttick

; middle part.
tpos = tposs[*,1]
yrange = lin_yrange
ylog = 0
ytickn = lin_ytickn
yspace = (yticklen ge 0)? 0: yticklen
tx = tpos[0]-xchsz*0.15+yspace
ty = tpos[1]-ychsz*0.25
xyouts, tx,ty,/normal,alignment=1, lin_ytickn[0]
ty = tpos[3]-ychsz*0.25
xyouts, tx,ty,/normal,alignment=1, lin_ytickn[lin_yticks]

ytickn[[0,lin_yticks]] = ' '
store_data, tvar, limits=lim
store_data, tvar, limits={$
    yrange: yrange, $
    ytitle: '('+unit+')', $
    labels: '', $
    xticklen: xticklen, $
    ytickv: lin_ytickv, $
    yticks: lin_yticks, $
    yminor: lin_yminor, $
    ytickname: ytickn, $
    yticklen: yticklen, $
    ylog: ylog}
tplot, tvar, trange=full_time_range, position=tpos, /noerase, /nouttick
tx = tpos[2]+xchsz*0.5
ty = (tpos[1]+tpos[3])*0.5
msg = strupcase(mission)
if mission ne 'polar' then msg += '-'+strupcase(probe)
msg += '!CUVW E!Du'
xyouts, tx,ty,/normal, msg
plots, tpos[[0,2]],(tpos[1]+tpos[3])*0.5+[0,0],/normal, color=white, linestyle=1


;---Secend panel.
tpos = poss[*,1]
tx = xchsz*0.5
ty = tpos[3]-ychsz*0.5
xyouts, tx,ty,/normal, 'b. Morlet'
options, amp_var, 'yticklen', yticklen
tplot, amp_var, trange=full_time_range, position=tpos, /noerase
get_data, amp_var, limits=lim
plot, full_time_range, lim.yrange, /nodata, /noerase, position=tpos, $
    xstyle=5, xlog=0, ystyle=5, ylog=1
foreach period, const_periods do begin
    tx = tpos[0]+xchsz*0.4
    tmp = convert_coord(full_time_range[0],period, /data, /to_normal)
    ty = tmp[1]-ychsz*0.2
    msg = string(round(spin_rate/period),format='(I0)')
    if msg eq '1' then msg = 'f' else msg += 'f'
    xyouts, tx,ty,/normal, alignment=0.5, msg, color=white
endforeach


if keyword_set(test) then stop
sgclose

end