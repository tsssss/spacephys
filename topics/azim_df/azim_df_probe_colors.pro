;+
; Return 6 colors for differentiate probes.
;-
function azim_df_probe_colors, nprobe

    ;probe_colors = sgcolor(['firebrick','orange','sea_green','deep_sky_blue','medium_blue','dark_violet'])
    if n_elements(nprobe) eq 0 then nprobe = 6
    probe_colors = smkarthm(90,240,nprobe, 'n')
    probe_color_ct = 64
    foreach color, probe_colors, ii do probe_colors[ii] = sgcolor(color,ct=probe_color_ct)

    return, probe_colors

end
