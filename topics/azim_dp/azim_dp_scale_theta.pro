;+
; Scale theta in terms of MLT, following the Gauss fit over all dipolarizations.
;
; theta. The theta data, a number or an array.
; mlt. The MLT data, a number of an array.
; width=. The scale width in MLT.
;-

function azim_dp_scale_theta, theta, mlt, width=width
    if n_elements(width) eq 0 then width = 4.
    z = (abs(mlt)<6)/width[0]
    return, theta*exp(z^2*0.5)
end
