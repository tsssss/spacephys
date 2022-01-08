;+
; read info in the log file, use the info to calc Poynting flux, and to 
; plot other data. Save plots into pdf files, save data into tplot file.
;-
pro cusp_plot_all_overview, eventid
    
    ; read info out of the log file.
    rootdir = shomedir()+'/Google Drive/works'
    if file_test(rootdir) eq 0 then rootdir = sdiskdir('Research')
    
    logfile = rootdir+'/works/cusp/cusp_list_of_conjun_9_10_all.log'
    info = cusp_read_conjun_list(logfile, event = eventid)
    if size(info,/type) ne 8 then begin
        print, 'no id found ...'
        return
    endif
    
    ; directory to save output data and pdfs.
    orootdir = shomedir()+'/cusp/overview'
    if ~file_test(orootdir,/directory) then file_mkdir, orootdir    
    
    
; **** basic info for both loading data and plotting.
    potr = info.polar.plot_time
    fatr = info.fast.plot_time
    potrcusp = info.polar.cusp_time
    fatrcusp = info.fast.cusp_time
    if potr[1] lt potr[0] then potr[1]+= 86400d
    if fatr[1] lt fatr[0] then fatr[1]+= 86400d
    if potrcusp[1] lt potrcusp[0] then potrcusp[1]+= 86400d
    if fatrcusp[1] lt fatrcusp[0] then fatrcusp[1]+= 86400d
    id = info.id


; **** plot settings and generate plots to disk.
    if keyword_set(no_plot) then return
    options, ['fa_','po_']+'dis', 'ytitle', 'R'
    rgb = [6,4,2] & red = 6
    charsz = 1
    !p.font = 1
    tplot_options, 'ygap', 0.25
    tplot_options, 'ynozero', 1
    tplot_options, 'version', 2
    tplot_options, 'num_lab_min', 8
    tplot_options, 'labflag', 1
    tplot_options, 'charsize', charsz*0.9
    tplot_options, 'xcharsize', charsz*0.7
    tplot_options, 'ycharsize', charsz*0.8
    tplot_options, 'zcharsize', charsz*0.6
    time_stamp, /off


    alltr = []
    ; use cusp time, expand 2.5 cusp time on both sides.
    ; polar ilat, mlt.
    j = 0
    dt = 600        ; 10 min.
    get_data, 'po_dis', t0, tmp
    if min(interpol(tmp, t0, potrcusp)) le 2.5 then dt = 120
    tmp = potrcusp
    ttr = 0.5*(tmp[1]+tmp[0])+[-1,1]*(tmp[1]-tmp[0])*2.5
    ttr = ttr-(ttr mod dt) & ttr[1]+= dt
    alltr = [min([ttr,alltr]),max([ttr,alltr])]
    ; fast ilat, mlt.
    j = 1
    dt = 60        ; 1 min.
    tmp = fatrcusp
    ttr = 0.5*(tmp[1]+tmp[0])+[-1,1]*(tmp[1]-tmp[0])*2.5
    ttr = ttr-(ttr mod dt) & ttr[1]+= dt
    alltr = [min([ttr,alltr]),max([ttr,alltr])]

    ; plot 7: omni, polar, fast overview.
    ofn = orootdir+'/'+id+'_overview.eps'
    sgpsopen, ofn, xsize = 6, ysize = 8, /inch
    plot_polar_fast_summary, stoepoch(alltr,'unix'), /no_delete
    sgpsclose, /pdf

end


;eventids = ['1998_0925_06','1998_0925_05']
eventids = [$
    '1997_0908_16', '1997_0914_13', '1997_1001_12', '1997_1009_15', '1997_1012_13', $
    '1997_1026_15', '1998_0905_21', '1998_0908_02', '1998_0908_20', '1998_0912_14', $
    '1998_0914_19', '1998_0918_13', '1998_0919_23', '1998_0919_24', '1998_0922_04', $
    '1998_0922_23', '1998_0923_17', '1998_0925_05', '1998_0925_06', '1998_0928_21', $
    '1998_0928_23', '1998_0929_16', '1998_1001_02', '1998_1001_03', '1998_1001_21', $
    '1998_1002_16', '1998_1004_20', '1998_1007_20', '1998_1008_14', '1998_1008_15', $
    '1998_1009_08', '1998_1009_09', '1998_1016_18', '1998_1017_12', '1998_1018_05', $
    '1998_1019_18', '1998_1020_10', '1998_1020_11', '1998_1021_05', '1998_1022_18', $
    '1998_1023_11', '1998_1024_23', '1998_1026_10', '1998_1029_10', '1998_1029_11', $
    '1998_1030_03', '1998_1030_05', '1998_1031_14', '1998_1031_15', '1999_0902_13', $
    '1999_0904_17', '1999_0905_10', '1999_0905_12', '1999_0908_12', '1999_0910_16', $
    '1999_0911_12', '1999_0912_06', '1999_0914_12', '1999_0917_09', '1999_0917_11', $
    '1999_0920_09', '1999_0920_11', '1999_0922_15', '1999_0923_09', '1999_0929_10', $
    '1999_0930_04', '1999_1001_15', '1999_1002_08', '1999_1008_08', '1999_1009_20', $
    '1999_1010_15', '1999_1011_08', '1999_1020_06', '1999_1020_08', '1999_1022_14', $
    '1999_1028_15', '1999_1101_08']

foreach eventid, eventids do $
    cusp_plot_all_overview, eventid
end
