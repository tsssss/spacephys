;+
; Type:
; 	function.
; 
; Name:
; 	tplot_tr.
; 
; Purpose:
; 	Return given variable's time range in unix time.
; 
; Parameters:
; 	vname, in, type = string, required.
;		Variable name.
; 
; Keywords:
; 	full, in, type = boolean, optional.
;		Set to return full time range.
; 
; Return:
; 	return, out, type = dblarr[2]/-1.
;		Time range in unix time.
;		-1 if no variable is found.
; 
; Example:
; 	tr = tplot_tr('AE', /full).
; 
; Notes:
; 	none.
; 
; Dependence:
; 	tdas.
; 
; Author:
; 	Sheng Tian.
; 
; History:
; 	2013-06-19, Sheng Tian, create.
;-

function tplot_tr, vname, full = full
    
    @tplot_com.pro
    
    if n_elements(data_quants) eq 0 then return, -1
    
    vnames = data_quants.name
    idx = where(vnames eq vname)
    if idx[0] eq -1 then begin
        message, 'can not find '+vname+' ...', /continue
        return, -1
    endif
    
    if keyword_set(full) then begin
        get_data, vname, data = tmp
        return, [tmp.x[0],tmp.x[-1]]
    endif
    
    return, data_quants[idx].trange

end
