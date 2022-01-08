
ids = cusp_id('1999_09')
nid = n_elements(ids)

rootdir = shomedir()+'/Google Drive/works'
datdir = rootdir+'/data/cusp'
outdir = shomedir()+'/deng'
if file_test(outdir,/directory) eq 0 then file_mkdir, outdir

for i = 0, nid-1 do begin
    ; read data.
    tfn = datdir+'/'+ids[i]+'_all_data.tplot'
    
    ; extract parallel poynting flux.
    tplot_restore, filename = tfn
    
    vars = ['pf_fac_mat_para_map','ilat','mlt','dis']
    vars = ['po_'+vars,'fa_'+vars]
    
    ofn = outdir+'/deng_conjun_pflux_'+ids[i]+'.tplot'
    tplot_save, vars, filename = ofn
    
endfor





end