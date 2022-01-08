;+
; Type: <+++>.
; Purpose: <+++>.
; Parameters: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Keywords: <+++>.
;   <+varname+>, <+in/out+>, <+datatype+>, <+req/opt+>. <+++>.
; Return: <+++>.
; Notes: <+++>.
; Dependence: <+++>.
; History:
;   <+yyyy-mm-dd+>, Sheng Tian, create.
;-

pro cusp_aux_modify_all_data

    rootdir = shomedir()+'/Google Drive/works/data/cusp'

    fns = rootdir+'/*_all_data.tplot'
    fns = file_search(fns, count = nfn)
    if nfn eq 0 then message, 'no data ...'

    for i = 0, nfn-1 do begin
        store_data, '*', /delete

        print, fns[i]
        tplot_restore, filename = fns[i]
        get_data, 'scidat', tmp, info

        ; do something.        

        tplot_save, filename = fns[i]
    endfor
end
