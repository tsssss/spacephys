;+
; Move certain type of plot from each event's directory to a directory of
; this type plot for all events.
;-


pro cusp_move_files, ids, indir = idir, outdir = odir, newdir = filetype1, file = filetype2


; create current type's directory.
dir = odir+'/'+filetype1
if ~file_test(dir, /directory) then file_mkdir, dir

if n_elements(ids) ne 0 then fns = idir+'/'+ids else fns = idir+'/*'
fns = file_search(fns)

foreach tfn, fns do begin
    id = file_basename(tfn)
    bfn = id+'_'+filetype2+'.pdf'    ; base filename.
    if ~file_test(tfn, /directory) then continue
    id = strmid(tfn, strlen(idir)+1)
    ifn = tfn+'/'+bfn
    if file_test(ifn) then begin
        ofn = dir+'/'+bfn
        if file_test(ofn) then file_delete, ofn
        file_copy, ifn, ofn
    endif
endforeach

end


; the root directory of event based tree.
idir = sdiskdir('GoogleDrive')+'/My Drive/works/works/cusp/cusp list conjun'
idir = sdiskdir('Research')+'/works/cusp/plots_using_mat'
;idir = shomedir()+'/cusp'

; the root directory of type based tree.
odir = shomedir()+'/cusp/survey'

; the type name of the directory to move files to.
filetype1 = 'survey'

; the type name in original file name.
filetype2 = 'field_and_poynt_freq_band'

filetype1 = 'overview'
filetype2 = 'overview'
;
;filetype1 = 'poynt_bands'
;filetype2 = 'field_and_poynt_freq_band'
;
;filetype1 = 'fast_bands'
;filetype2 = 'mat_spec_fast'
;
;;filetype1 = 'polar_bands'
;;filetype2 = 'mat_spec_polar'
;
;filetype1 = 'hydra'
;filetype2 = 'hydra'
;
;filetype1 = 'esa'
;filetype2 = 'esa'
;
;filetype1 = 'ebcorr'
;filetype2 = 'polar_e56dot0_dbcorr'
;
;filetype1 = 'spc'
;filetype2 = 'polar_spc2fac'

ids = cusp_id_new('south_imf')


cusp_move_files, ids, indir = idir, outdir = odir, newdir = filetype1, file = filetype2
end