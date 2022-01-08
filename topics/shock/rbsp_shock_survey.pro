;+
; Type: procedure.
; Purpose: Survey plot to look at shock.
; Parameters:
;   trs, in, dblarr[n,2], req. Time ranges in ut.
; Keywords: none.
; Notes: none.
; Dependence: tdas,slib.
; History:
;   2014-02-21, Sheng Tian, create.
;-
pro rbsp_shock_survey, trs
    compile_opt idl2

    ; prepare.
    probes = ['a','b']
    nprobe = n_elements(probes)
    ntr = n_elements(trs)/2
    tpad = 1200D    ; 20 min.
    re = 6374D & re1 = 1D/re
    rbsp_efw_init
    istp_init
    pre1 = 'OMNI_HRO_1min_'

    ; loop for each time range.
    for i = 0, ntr-1 do begin
        ttr = [trs[i,0]-tpad,trs[i,1]+tpad]
        timespan, ttr, ttr[1]-ttr[0], /second
        if ttr[1]-ttr[0] gt 86400L*5 then message, 'more than 5 days ...'
        rbsp_load_spice_kernels
        
        ; load omni bz.
        omni_hro_load, trange = ttr
        stplot_renew, pre1+'BZ_GSM', newname = 'omni_bz_gsm'
        stplot_renew, pre1+'SYM_H', newname = 'omni_symh'
        store_data, pre1+'*', /delete

        ; loop for probes.
        for j = 0, nprobe-1 do begin
            pre0 = 'rbsp'+probes[j]+'_'
            ; load pos.
            rbsp_load_spice_state, probe = probes[j], coord = 'gse', $
                /no_spice_load
            get_data, pre0+'state_pos_gse', data = tmp
            tmp.y *= re1
            store_data, pre0+'pos_gse', data = tmp
            store_data, pre0+'state_*', /delete
            ; load de.
            rbsp_load_efw_waveform_l2, probe = probes[j], trange = ttr
            get_data, pre0+'efw_e-spinfit-mgse_e12_spinfit_mgse', data = tmp
            tmp.y[*,0] = 0 & store_data, pre0+'de_mgse', data = tmp
            rbsp_mgse2gse, pre0+'de_mgse', newname = pre0+'de_gse', $
                probe = probes[j], /no_spice_load
            store_data, pre0+'efw_*', /delete
            ; load b.
            rbsp_load_emfisis, probe = probes[j], trange = ttr, coord = 'gse'
            stplot_renew, pre0+'emfisis_l3_4sec_gse_Mag', newname = pre0+'b_gse'
            store_data, pre0+'emfisis_*', /delete
            get_data, pre0+'b_gse', tmp, db
            for k = 0, 2 do $
                db[*,k] = db[*,k]-smooth(db[*,k],30*60/sdatarate(tmp))
            store_data, pre0+'db_gse', tmp, db
            ylim, pre0+'db_gse', -10, 10, 0
            ; prepare vars.
            vars = pre0+['pos_gse','b_gse','de_gse','db_gse']
            options, vars, 'colors', [6,4,2]
            options, vars, 'labels', ['x','y','z']
            vars = [vars,'omni_'+['bz_gsm','symh']]
            tplot, vars, trange = trs[i,*]
            stop
        endfor
    endfor
end

trs = time_double([['2012-09-30/00:00'],['2012-10-01/00:00']])
trs = time_double([['2012-10-08/00:00'],['2012-10-09/00:00']])
trs = time_double([['2012-10-12/00:00'],['2012-10-13/00:00']])
trs = time_double([['2012-11-01/00:00'],['2012-11-02/00:00']])
trs = time_double([['2012-11-14/00:00'],['2012-11-15/00:00']])
trs = time_double([['2012-11-19/00:00'],['2012-11-20/00:00']])
trs = time_double([['2012-11-23/00:00'],['2012-11-24/00:00']])
trs = time_double([['2013-01-16/00:00'],['2013-01-17/00:00']])
trs = time_double([['2013-01-23/00:00'],['2013-01-24/00:00']])
trs = time_double([['2013-03-17/00:00'],['2013-03-18/00:00']])
trs = time_double([['2013-03-20/00:00'],['2013-03-21/00:00']])
trs = time_double([['2013-03-29/00:00'],['2013-03-30/00:00']])
trs = time_double([['2013-04-24/00:00'],['2013-04-25/00:00']])
trs = time_double([['2013-04-30/00:00'],['2013-05-01/00:00']])
trs = time_double([['2013-05-31/00:00'],['2013-06-01/00:00']])
trs = time_double([['2013-05-25/00:00'],['2013-05-26/00:00']])
trs = time_double([['2013-06-01/00:00'],['2013-06-02/00:00']])
trs = time_double([['2013-06-06/00:00'],['2013-06-07/00:00']])
trs = time_double([['2013-06-27/00:00'],['2013-06-28/00:00']])
trs = time_double([['2013-07-06/00:00'],['2013-07-07/00:00']])
trs = time_double([['2013-07-09/00:00'],['2013-07-11/00:00']])
trs = time_double([['2013-07-15/00:00'],['2013-07-16/00:00']])
trs = time_double([['2013-08-04/00:00'],['2013-08-05/00:00']])
trs = time_double([['2013-08-13/00:00'],['2013-08-14/00:00']])
trs = time_double([['2013-08-27/00:00'],['2013-08-28/00:00']])
trs = time_double([['2013-08-30/00:00'],['2013-08-31/00:00']])
trs = time_double([['2013-10-02/00:00'],['2013-10-03/00:00']])
trs = time_double([['2013-10-08/00:00'],['2013-10-09/00:00']])
trs = time_double([['2013-10-14/00:00'],['2013-10-15/00:00']])
trs = time_double([['2013-10-30/00:00'],['2013-10-31/00:00']])
trs = time_double([['2013-11-07/00:00'],['2013-11-08/00:00']])
trs = time_double([['2013-11-09/00:00'],['2013-11-10/00:00']])
trs = time_double([['2013-11-11/00:00'],['2013-11-12/00:00']])
trs = time_double([['2013-12-07/00:00'],['2013-12-08/00:00']])
trs = time_double([['2013-12-08/00:00'],['2013-12-09/00:00']])
trs = time_double([['2013-12-13/00:00'],['2013-12-14/00:00']])
trs = time_double([['2013-12-25/00:00'],['2013-12-26/00:00']])
trs = time_double([['2014-02-07/00:00'],['2014-02-08/00:00']])
trs = time_double([['2014-02-08/00:00'],['2014-02-09/00:00']])
trs = time_double([['2014-02-15/00:00'],['2014-02-16/00:00']])
trs = time_double([['2014-02-18/00:00'],['2014-02-19/00:00']])
trs = time_double([['2014-02-20/00:00'],['2014-02-21/00:00']])
trs = time_double([['2014-02-27/00:00'],['2014-02-28/00:00']])

trs = time_double([['2013-06-06/00:00'],['2013-06-07/00:00']])

rbsp_shock_survey, trs
end
