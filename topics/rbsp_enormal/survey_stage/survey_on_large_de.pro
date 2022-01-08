;+
; plot key quantities on screen.
; determine event id, start/end times, dst, ae, has hope spec, has jump in beta, density.
;-

pro survey_on_large_de, utr, tprobe, logfn = ologfn

    tinfo = {id:'', $           ; event id YYYY_MMDD_hhmm.
        tsta:0d, tend:0d, $     ; star/end time of the event.
        dst:0d, ae:0d, $        ; min Dst, max AE.
        ebipolar:-1d, $                         ; 0/1/-1 for no/yes/no value.
        ezmgse:0d, ezdot0mgse:0d, $             ; max abs value (with sign), defined when ebipolar is 0.
        jumpinbeta:-1d, jumpindensity:-1d, $    ; 0/1/-1 for no/yes/no value.
        bzincrease:-1d, bincrease:-1, $         ; 0/1/-1 for no/yes/no value.
        fptmlat:0d, fptmlt:0d, mlt:0d}          ; location info.


    pre0 = 'rbsp'+tprobe+'_'
    deg = 180d/!dpi
    rad = !dpi/180
    re = 6378d & re1 = 1d/re
    rgb = [6,4,2]
    dt0 = 86400d
    ratio0 = 1.5        ; threshold for dE asymmetry, bipolar if smaller than ratio0.
    utr0 = utr-(utr mod dt0)+[0,dt0]
    timespan, utr0[0], dt0, /second


    tplot_options, 'constant', 0
    tplot_options, 'num_lab_min', 10
    tplot_options, 'labflag', -1
    tplot_options, 'ynozero', 1
    tplot_options, 'ymargin', [5,5]
    tplot_options, 'xmargin', [25,15]


; **** load omni data.
; dst, ae.
    omni = sread_omni(utr0)
    if size(omni,/type) eq 8 then begin
        uts = sfmepoch(omni.epoch,'unix')
        store_data, 'dst', uts, omni.symh, limits = {ytitle:'Dst (nT)'}
        store_data, 'ae', uts, omni.ae, limits = {ytitle:'AE (nT)'}
    endif

; **** load hope data.
; rbspx_[h,o,e]_en, rbspx_[h,o,e]_pa, rbspx_[n,t].
    vars = ['PITCH_ANGLE',$
            'Epoch_Ele','HOPE_ENERGY_Ele',$
            'Epoch_Ion','HOPE_ENERGY_Ion',$
            'FEDU','FPDU','FODU']
    hopel2 = sread_rbsp_hope_l3(utr0, probes = tprobe, vars = vars)
    if size(hopel2,/type) eq 8 then begin
        uts = sfmepoch(hopel2.epoch_ion,'unix')
        store_data, pre0+'h_en', uts, total(hopel2.fpdu,3,/nan), hopel2.hope_energy_ion
        store_data, pre0+'h_pa', uts, total(hopel2.fpdu,2,/nan), hopel2.pitch_angle
        store_data, pre0+'o_en', uts, total(hopel2.fodu,3,/nan), hopel2.hope_energy_ion
        store_data, pre0+'o_pa', uts, total(hopel2.fodu,2,/nan), hopel2.pitch_angle
        uts = sfmepoch(hopel2.epoch_ele,'unix')
        store_data, pre0+'e_en', uts, total(hopel2.fedu,3,/nan), hopel2.hope_energy_ele
        store_data, pre0+'e_pa', uts, total(hopel2.fedu,2,/nan), hopel2.pitch_angle

        vars = pre0+['h_en','h_pa','o_en','o_pa','e_en','e_pa']
        options, vars, 'spec', 1
        options, vars, 'no_interp', 1
        options, vars, 'zlog', 1
        options, vars, 'ztitle', '(s!E-1!Ncm!E-2!Nsr!E-1!NkeV!E-1!N)'

        vars = pre0+['h_en','o_en','e_en']
        options, vars, 'ylog', 1
        options, vars, 'yrange', [10,4e4]
        
        vars = pre0+['e_en','e_pa']
        options, vars, 'zrange', [1e4,1e10]
        vars = pre0+['h_en','h_pa']
        options, vars, 'zrange', [1e4,1e7]
        vars = pre0+['o_en','o_pa']
        options, vars, 'zrange', [1e4,1e6]

        vars = pre0+['h_pa','o_pa','e_pa']
        options, vars, 'yrange', [0,180]
        
        tvar = pre0+'e_en'
        options, tvar, 'ytitle', 'e- energy!C(eV)'
        tvar = pre0+'h_en'
        options, tvar, 'ytitle', 'H+ energy!C(eV)'
        tvar = pre0+'o_en'
        options, tvar, 'ytitle', 'O+ energy!C(eV)'

        tvar = pre0+'e_pa'
        options, tvar, 'ytitle', 'e- pitch!C(deg)'
        tvar = pre0+'h_pa'
        options, tvar, 'ytitle', 'H+ pitch!C(deg)'
        tvar = pre0+'o_pa'
        options, tvar, 'ytitle', 'O+ pitch!C(deg)'
    endif
    
    hopemom = sread_rbsp_hope_l3(utr0, probes = tprobe, type = 'mom')
    if size(hopemom,/type) eq 8 then begin
        uts = sfmepoch(hopemom.epoch_ele,'unix')
        store_data, pre0+'n', uts, hopemom.dens_e_200, $
            limits = {ytitle:'N!Ie!N!C(cm!E-3!N)', ylog:1, constant:1, labels:'Ne'}
        store_data, pre0+'t', uts, [[hopemom.tperp_e_200],[hopemom.tpar_e_200]], $
            limits = {ytitle:'T!I!N!C(eV)', ylog:1, colors:[6,0], labels:['Tperp','Tpara'], constant:1000}
        tvar = pre0+'t'
        get_data, tvar, uts, dat
        idx = where(dat eq 1e20, cnt)
        if cnt ne 0 then dat[idx] = !values.d_nan
        store_data, tvar, uts, dat
    endif


; **** load E field. l3, mgse. And pos_gse, mlt/mlat/lshell.
; rbspx_de_[mgse,svy], rbspx_[mlt,lshell,ilat,mlat].
    efwl3 = sread_rbsp_efw_l3(utr, probes = tprobe)
    uts = sfmepoch(efwl3.epoch,'unix',/epoch16)

    store_data, pre0+'de_mgse', uts, efwl3.efield_inertial_frame_mgse, $
        limits = {ytitle:'dE!C(mV/m)', colors:rgb, labels:'MGSE '+['x','y','z']}
    store_data, pre0+'mlt', uts, (efwl3.mlt_lshell_mlat)[*,0], $
        limits = {ytitle:'MLT (hr)'}
    store_data, pre0+'lshell', uts, (efwl3.mlt_lshell_mlat)[*,1], $
        limits = {ytitle:'L (Re)'}
    store_data, pre0+'ilat', uts, acos(1/sqrt((efwl3.mlt_lshell_mlat)[*,1]))*deg, $
        limits = {ytitle:'ILat (deg)'}
    store_data, pre0+'mlat', uts, (efwl3.mlt_lshell_mlat)[*,2], $
        limits = {ytitle:'MLat (deg)'}
    store_data, pre0+'pos_gse', uts, efwl3.pos_gse, $
        limits = {colors:rgb, labels:['x','y','z']}

    
    vars = ['esvy']
    rbsp_load_efw_waveform, probe = tprobe, trange = utr0, datatype = vars, type = 'calibrated'
    tvar = pre0+'efw_esvy'
    get_data, tvar, tmp, dat
    dat[*,2] = 0
    store_data, pre0+'de_svy', tmp, dat, $
        limits = {ytitle:'dE survey!C(mV/m)', colors:rgb, labels:'UVW '+['x','y','z']}
    store_data, pre0+'efw_esvy*', /delete

    
; **** load B field. l3, gse, emfisis, 1sec.
; rbspx_[b,bz]_gse, rbspx_b.
    emfisis = sread_rbsp_emfisis(utr0, probes = tprobe)
    uts = sfmepoch(emfisis.epoch,'unix',/tt2000)
    store_data, pre0+'b_gse', uts, emfisis.mag, $
        limits = {ytitle:'B!C(nT)',colors:rgb, labels:'GSE '+['x','y','z']}
    store_data, pre0+'b', uts, sqrt(total(emfisis.mag^2,2)), $
        limits = {ytitle:'B total!C(nT)', labels:'B total'}
    store_data, pre0+'bz_gse', uts, emfisis.mag[*,2], $
        limits = {ytitle:'Bz!C(nT)', labels:'Bz GSE'}


; **** map footprint.
; rbspx_fpt_[lat,lon].
    ; prepare mapping.
    r0 = 1+110*re1  ; 110 km altitude.
    model = 't89'
    dir = -1        ; always north hem, b/c conjugate to thm_asi.
    r0 = 1+110*re1  ; 110 km altitude.
    sgeopack_par, utr, model, /delete  ; get tplot var <model>_par.
    t89 = 0 & t96 = 0 & t01 = 0
    case model of
        't89': t89 = 1
        't96': t96 = 1
        't01': t01 = 1
    endcase
    ; map pos to fpt, convert to mag.
    get_data, pre0+'pos_gse', data = tmp
    uts = tmp.x & ets = 1000D*uts+62167219200000D
    pos0 = tmp.y*re1 & pos1 = pos0      ; in re.
    ; interpolate par.
    if model ne 't89' then begin
        get_data, model+'_par', data = tmp
        pars = sinterpol(tmp.y, tmp.x, uts)
    endif
    ; loop for each time.
    for j = 0, n_elements(uts)-1 do begin
        ; set geopack.
        geopack_epoch, ets[j], yr, mo, dy, hr, mi, sc, msc, /breakdown_epoch
        geopack_recalc, yr, mo, dy, hr, mi, sc+msc*0.001D, /date
        ; pos in gse, which is the mapping coord.
        x0 = pos0[j,0] & y0 = pos0[j,1] & z0 = pos0[j,2]
        ; convert from gse to gsm.
        geopack_conv_coord, x0, y0, z0, /from_gse, $
            x1, y1, z1, /to_gsm
        dir = (z1 gt 0)? -1: 1
        if model ne 't89' then par = reform(pars[j,*]) else par = 2
        geopack_trace, x1, y1, z1, dir, par, xf, yf, zf, $
            epoch = ets[j], /refine, /ionosphere, $
            t89 = t89, t96 = t96, t01 = t01
        ; convert from gse to mag.
        geopack_conv_coord, xf, yf, zf, /from_gsm, $
            x1, y1, z1, /to_mag
        pos1[j,*] = [x1,y1,z1]
    endfor
    mlat = asin(pos1[*,2]*(1/r0))*deg
    mlon = atan(pos1[*,1],pos1[*,0])*deg
    store_data, pre0+'fpt_mag', data = {x:uts, y:pos1}
    store_data, pre0+'fpt_mlat', data = {x:uts, y:mlat}, $
        limits = {ytitle:'MLat/fpt (deg)'}
    store_data, pre0+'fpt_mlt', uts, slon2lt(mlon, stoepoch(uts,'unix'), /mag, /deg)/15, $
        limits = {ytitle:'MLT/fpt (hr)'}

    
; **** other variables.
; rbspx_beta, rbspx_[b,de_dot0]_[gse,mgse].
    ; plasma beta. 2mu0*nkT/B^2
    tmp = 2*4*!dpi*1e-7*1e6*1.6e-19*1e18    ; 
    get_data, pre0+'n', t0, dat & dat*= tmp
    get_data, pre0+'t', t0, tmp & dat*= tmp[*,0]
    get_data, pre0+'b', t1, tmp & tmp = interpol(tmp, t1, t0) & dat/= tmp^2
    store_data, pre0+'beta', t0, dat, $
        limits = {ytitle:'Beta', labels:'Beta'}

;    get_data, pre0+'de_gse', t0, de_mgse
;    de_mgse[*,0] = 0
;    store_data, pre0+'de_mgse', t0, de_mgse
;    rbsp_mgse2gse, pre0+'de_mgse', newname = pre0+'de_gse'
;    get_data, pre0+'de_gse', t0, de_gse
;    get_data, pre0+'b_gse', tmp, b_gse
;    b_gse = sinterpol(b_gse, tmp, t0)
;    store_data, pre0+'b_gse', t0, b_gse
;    rbsp_mgse2gse, pre0+'b_gse', newname = pre0+'b_mgse', /inverse
;    get_data, pre0+'b_mgse', tmp, b_mgse
;    de_mgse[*,0] = -(b_mgse[*,1]*de_mgse[*,1]+b_mgse[*,2]*de_mgse[*,2])/b_mgse[*,0]
;    ratio = abs(b_mgse[*,0])/snorm(b_mgse)
;    store_data, pre0+'bx_ratio', t0, ratio
;    idx = where(ratio le minratio, cnt)
;    if cnt ne 0 then de_mgse[idx,0] = !values.d_nan
;    store_data, pre0+'de_dot0_mgse', t0, de_mgse
;    rbsp_mgse2gse, pre0+'de_dot0_mgse', newname = pre0+'de_dot0_gse'


; **** start survey on screen.
    vars = pre0+['de_mgse','de_svy','b','bz_gse','e_en','h_en','n','t','beta']
    labs = [pre0+['mlt','lshell','ilat','fpt_mlat','fpt_mlt'],'dst','ae']
    titl = 'RBSP-'+strupcase(tprobe)+' survey plot, '+ $
        time_string(utr[0])+' to '+time_string(utr[1],tformat='hh:mm')
    
    sgopen, 0, xsize = 800, ysize = 900
    device, decomposed = 0
    loadct2, 43
    tplot, vars, var_label = labs, trange = utr, title = titl
    
    msg = ''
    read, msg, prompt = 'Zoom in on this plot? (y/n)'
    while msg eq 'y' do begin
        ctime, tutr, /exact
        while n_elements(tutr) ne 2 do ctime, tutr, /exact
        tplot, vars, var_label = labs, trange = tutr, title = titl
        read, msg, prompt = 'Zoom out? (y/n)'
        if msg eq 'y' then begin        ; change to a different time.
            tplot, vars, var_label = labs, trange = utr, title = titl
            read, msg, prompt = 'Zoom in on this plot? (y/n)'
            ctime, tutr, /exact
        endif else begin                ; work on the current time.
            read, msg, prompt = 'Set output: (d/b/m), d:dipolarization, b:boundary crossing, m:manual set. '
            while msg ne 'd' and msg ne 'b' and msg ne 'm' do $
                read, msg, prompt = 'Set output: (d/b/m), d:dipolarization, b:boundary crossing, m:manual set. '
            case msg of
                'd': begin
                    tinfo.bzincrease = 1
                    tinfo.jumpinbeta = 1
                    tinfo.jumpindensity = 1
                    end
                'b': begin
                    tinfo.bzincrease = 0
                    tinfo.jumpinbeta = 1
                    tinfo.jumpindensity = 1
                    end
                'm': begin
                    tmp = 0
                    read, tmp, prompt = 'Bz increase? (0/1/-1)'
                    tinfo.bzincrease = tmp 
                    read, tmp, prompt = 'Jump in beta? (0/1/-1)'
                    tinfo.jumpinbeta = tmp 
                    read, tmp, prompt = 'Jump in density? (0/1/-1)'
                    tinfo.jumpindensity = tmp 
                    end
            endcase
            
            ; prepare info.
            tinfo.id = time_string(tutr[0],tformat='YYYY_MMDD_hhmm')
            tinfo.tsta = tutr[0]-(tutr[0] mod 60)
            tinfo.tend = tutr[1]-(tutr[1] mod 60)+60

            get_data, 'dst', t0, dat
            idx = where(t0 ge tutr[0] and t0 le tutr[1], cnt)
            tinfo.dst = (cnt eq 0)? !values.d_nan: min(dat[idx])

            get_data, 'ae', t0, dat
            idx = where(t0 ge tutr[0] and t0 le tutr[1], cnt)
            tinfo.ae = (cnt eq 0)? !values.d_nan: max(dat[idx])

            get_data, pre0+'de_mgse', t0, dat
            idx = where(t0 ge tutr[0] and t0 le tutr[1], cnt)
            if cnt eq 0 then begin
                tinfo.ezmgse = !values.d_nan
                tinfo.ezdot0mgse = !values.d_nan
                tinfo.ebipolar = -1
            endif else begin
                maxde = max(dat[idx,2])
                minde = min(dat[idx,2])
                tmp = abs(maxde/minde) & if tmp lt 1 then tmp = 1d/tmp
                if tmp ge ratio0 then begin
                    tinfo.ezmgse = max([maxde,minde],/absolute)
                    tinfo.ezdot0mgse = !values.d_nan
                    tinfo.ebipolar = 0
                endif else begin
                    tinfo.ezmgse = !values.d_nan
                    tinfo.ezdot0mgse = !values.d_nan
                    tinfo.ebipolar = 1
                endelse
            endelse

            get_data, pre0+'fpt_mlat', t0, dat
            idx = where(t0 ge tutr[0] and t0 le tutr[1], cnt)
            tinfo.fptmlat = (cnt eq 0)? !values.d_nan: mean(dat[idx])

            get_data, pre0+'fpt_mlt', t0, dat
            idx = where(t0 ge tutr[0] and t0 le tutr[1], cnt)
            tinfo.fptmlt = (cnt eq 0)? !values.d_nan: mean(dat[idx])

            get_data, pre0+'mlt', t0, dat
            idx = where(t0 ge tutr[0] and t0 le tutr[1], cnt)
            tinfo.mlt = (cnt eq 0)? !values.d_nan: mean(dat[idx])

            tmp = 0
            read, tmp, prompt = 'B total increase? (0/1/-1)'
            tinfo.bincrease = tmp 

            
            ; output to log file.
            spc1 = '    '
            spc2 = '  '

            cmd = ''
            cmd+= tinfo.id+'_'+tprobe+spc1
            cmd+= time_string(tinfo.tsta,tformat='hh:mm')+spc1
            cmd+= time_string(tinfo.tend,tformat='hh:mm')+spc1
            cmd+= (finite(tinfo.dst,/nan)? 'nan': string(tinfo.dst,format='(I4)'))+spc1
            cmd+= (finite(tinfo.ae,/nan)? 'nan': string(tinfo.ae,format='(I4)'))+spc1
            case tinfo.ebipolar of
                0: tmp = ' n '
                1: tmp = ' y '
                -1: tmp = 'n/a'
            endcase
            cmd+= tmp+spc1
            cmd+= (finite(tinfo.ezmgse,/nan)? '  nan': string(tinfo.ezmgse,format='(F5.1)'))+spc2
            cmd+= (finite(tinfo.ezdot0mgse,/nan)? '  nan': string(tinfo.ezdot0mgse,format='(F5.1)'))+spc1

            cmd+= (finite(tinfo.fptmlat,/nan)? '  nan': string(tinfo.fptmlat,format='(F5.1)'))+spc2
            cmd+= (finite(tinfo.fptmlt,/nan)? '  nan': string(tinfo.fptmlt,format='(F5.1)'))+spc1
            cmd+= (finite(tinfo.mlt,/nan)? '  nan': string(tinfo.mlt,format='(F5.1)'))+spc1

            case tinfo.bzincrease of
                0: tmp = ' n '
                1: tmp = ' y '
                -1: tmp = 'n/a'
            endcase
            cmd+= tmp+spc2
            case tinfo.bincrease of
                0: tmp = ' n '
                1: tmp = ' y '
                -1: tmp = 'n/a'
            endcase
            cmd+= tmp+spc1

            case tinfo.jumpinbeta of
                0: tmp = ' n '
                1: tmp = ' y '
                -1: tmp = 'n/a'
            endcase
            cmd+= tmp+spc2
            case tinfo.jumpindensity of
                0: tmp = ' n '
                1: tmp = ' y '
                -1: tmp = 'n/a'
            endcase
            cmd+= tmp+spc1

            hd0 = ''
            hd0+= '    event id    '+spc1
            hd0+= 'start/end time'+spc1
            hd0+= 'Symh'+spc1
            hd0+= ' AE '+spc1
            hd0+= 'bi-'+spc1
            hd0+= 'E max (mV/m)'+spc1
            hd0+= 'FPT MLat/MLT'+spc1
            hd0+= ' MLT '+spc1
            hd0+= 'IncrFlag'+spc1
            hd0+= 'JumpFlag'+spc1

            hd1 = ''
            hd1+= 'YYYY_MMDD_hhmm'+spc1
            hd1+= 'hh:mm'+spc1
            hd1+= 'hh:mm'+spc1
            hd1+= '(nT)'+spc1
            hd1+= '(nT)'+spc1
            hd1+= 'plr'+spc1
            hd1+= ' dEz '+spc2
            hd1+= 'dEz_0'+spc1
            hd1+= '(deg)'+spc2
            hd1+= '(hrs)'+spc1
            hd1+= '(hrs)'+spc1
            hd1+= ' bz'+spc2
            hd1+= ' b '+spc1
            hd1+= 'beta'+spc2
            hd1+= 'n '+spc1


            hd2 = ''
            hd2+= '----------------'+spc1
            hd2+= '-----'+spc1
            hd2+= '-----'+spc1
            hd2+= '----'+spc1
            hd2+= '----'+spc1
            hd2+= '---'+spc1
            hd2+= '-----'+spc2
            hd2+= '-----'+spc1
            hd2+= '-----'+spc2
            hd2+= '-----'+spc1
            hd2+= '-----'+spc1
            hd2+= '---- ---'+spc1
            hd2+= '---- ---'+spc1

            header = [hd0,hd1,hd2]
            
            print, [header,cmd]

            if file_test(ologfn) eq 0 then begin
                file_mkdir, file_dirname(ologfn)
                stouch, ologfn
                addheader = 1
            endif else addheader = 0
            openw, lun, ologfn, /get_lun, /append
            if addheader then printf, lun, header
            printf, lun, cmd
            free_lun, lun
        endelse
        
        read, msg, prompt = 'Zoom out? (y/n)'
        if msg eq 'y' then begin        ; change to a different time.
            tplot, vars, var_label = labs, trange = utr, title = titl
            read, msg, prompt = 'Zoom in on this plot? (y/n)'
        endif
    endwhile 
    
end


ilogfn = shomedir()+'/Google Drive/works/works/rbsp_de/dipolarization/list_large_de.log'
ologfn = shomedir()+'/list_large_de_round3.log'

nheader = 3
headers = strarr(nheader)
nline = file_lines(ilogfn)-nheader
lines = strarr(nline)
openr, lun, ilogfn, /get_lun
readf, lun, headers
readf, lun, lines
free_lun, lun

for i = 106, nline-1 do begin
    tfn = file_basename(lines[i])
    tprobe = strlowcase(strmid(tfn,19,1))
    id = strmid(tfn,0,14)

    ut0 = time_double(id, tformat='YYYY_MMDD_hhmm')
    utr = ut0+[-1,1]*15*60

    survey_on_large_de, utr, tprobe, logfn = ologfn
endfor

end
