
; check RBSP s/c potential and density.


utr = time_double(['2012-11-14/04:41','2012-11-14/04:46'])
tprobe = 'a'


pre0 = 'rbsp'+tprobe+'_'

load = tnames(pre0+'vsc') eq ''
load = 1
if load then begin
    tmp = sread_rbsp_efw_l1(utr, probes = tprobe, type = ['vsvy'], vars = ['epoch','vsvy'])
    uts = sfmepoch(tmp.epoch,'unix')
    idx = where(uts ge utr[0] and uts le utr[1])
    uts = uts[idx]
    vsc = tmp.vsvy[idx,*]*0.01
    store_data, pre0+'v', uts, vsc, limits = $
        {ytitle:'V (V)', colors:[1,2,3,4,5,6], labels:'V'+['1','2','3','4','5','6']}
        
    vsc = mean(vsc[*,0:3], dimension = 2)
    store_data, pre0+'vsc', uts, vsc, limits = $
        {ytitle:'Vsc (V)', labels:'Vsc'}
endif

    
; convert Vsc to density, from Scott.
get_data, pre0+'vsc', uts, vsc
if tprobe eq 'a' then begin
    den =  3448.12*exp(vsc*2.75) +64.372*exp(vsc*0.474)
endif else begin
    den = 1818.5466*exp(vsc*2.4191940)+19.646707*exp(vsc*0.28274512)
endelse
den = 0.15*exp(vsc*0.045)   ; corresponds to 22 eV electrons.
store_data, pre0+'efw_ne', uts, den, limits = $
    {ytitle:'EFW N!De!C(cm!U-3!N)', ylog:1, constant:[0.1,1], labels:'EFW N!Dele'}
    

; load HOPE density.
hopemom = sread_rbsp_hope_l3(utr, probes = tprobe, type = 'mom')
uts = sfmepoch(hopemom.epoch_ele,'unix')
store_data, pre0+'hope_ne', uts, hopemom.dens_e_200, $
    limits = {ytitle:'HOPE N!De!N!C(cm!U-3!N)', ylog:1, constant:[0.1,1], labels:'HOPE N!Dele'}

store_data, pre0+'ne', data = pre0+['efw_ne','hope_ne'], limits = {colors:[6,0],labels:['','']}
    
tplot_options, 'labflag', -1
tplot_options, 'xmargin', [25,10]

ofn = shomedir()+'/test_rb_vsc_density.pdf'
ofn = 0
sgopen, ofn, xsize = 11, ysize = 8.5, /inch
device, decomposed = 0
loadct2, 43

vars = pre0+['v','vsc','efw_ne','hope_ne','ne']
tplot, vars, trange = utr

sgclose


end
