;+
; Preprocess E field in GSM: remove DC background.
;-

pro stplot_prep_efield, evar, newname=newname, $
    addto=varlist

    rgb = [6,4,2]
    xyz = ['x','y','z']
    
    pre0 = strmid(evar,0,strpos(evar,'_'))+'_'  ; '' if no _.
    get_data, evar, uts, egsm
    nrec = n_elements(uts)

    for i=0,2 do begin
        egsm[*,i] -= mean(egsm[*,i])
    endfor
    store_data, pre0+'de_gsm', uts, egsm, limits=$
        {ytitle:'(mV/m)', colors:rgb, labels:'GSM dE'+xyz, labflag:-1}
    
    myvars = pre0+['de_gsm']
    if keyword_set(newname) then myvars = newname
    if not keyword_set(varlist) then varlist=[]
    varlist = [varlist, myvars]

end
