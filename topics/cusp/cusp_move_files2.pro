
idir = shomedir()+'/cusp'
odir = shomedir()+'/Google Drive/works/works/cusp/cusp list conjun'
fns = file_search(idir+'/*', count = nfn)

thefns = ['hydra.pdf','ke_special.pdf','ke_special.svg']
thefns = ['efluxes','esa','field_and_poynt','field_and_poynt_freq_band'+['','_scaled'], $
    'mat_spec_'+['fast','polar'],'overview']+'.pdf'

for i = 0, nfn-1 do begin
    tfn = fns[i]
    tid = file_basename(tfn)
    for j = 0, n_elements(thefns)-1 do begin
        ifn = tfn+'/'+tid+'_'+thefns[j]
        ofn = odir+'/'+tid+'/'+tid+'_'+thefns[j]
;        ofn = odir+'/papco_backup/'+tid+'_'+thefns[j]
        if file_test(ifn) then file_copy, ifn, ofn, /overwrite
    endfor
endfor

end
