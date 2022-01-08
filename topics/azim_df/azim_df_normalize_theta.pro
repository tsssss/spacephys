;+
; Return color after normalizing theta (detrended tilt angle) according to MLT.
; Adopted from normalize_data.
;-

function azim_df_normalize_theta, zzs, mlt, smooth_width=smooth_width, zrange=range, ct=ct, reverse_ct=reverse_ct

    if n_elements(smooth_width) ne 0 then begin
        zzs -= smooth(zzs,smooth_width, /edge_truncate, /nan)
        mean = mean(zzs,/nan)
        zzs = zzs-mean
    endif
    if n_elements(mlt) ne 0 then zzs *= sqrt(abs(mlt)+1)

    linear_range = 1
    signs = sign(zzs)
    data = alog10(abs(zzs)>linear_range)*signs
    the_range = [-1,1]*alog10(max(abs(range)))
    index_color = bytscl(data, min=the_range[0], max=the_range[1])
    if keyword_set(reverse_ct) then index_color = 256-index_color

    if n_elements(ct) eq 0 then ct = 70
    loadct, ct, rgb_table=rgb0
    rgb0 = float(rgb0)
    true_colors = 256L*(256L*rgb0[index_color,2]+rgb0[index_color,1])+rgb0[index_color,0]

    return, true_colors
end
