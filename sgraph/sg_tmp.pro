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

pro sg_tmp

device, get_visual_name = vclass, get_visual_depth = vdepth
print, 'Current viaual class: '+vclass+', '+string(vdepth,format='(I2)'), '-bit'
if vdepth eq 24 then begin
    print, 'Please add the following command to IDL startup file: '
    print, '    device, true_color = 24'
    return
endif

end
