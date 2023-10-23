;+
; Plot RBSP-A and -B wave data.
;-

function fig_2015_0416_0800_wave_v01, plot_file, event_info=event_info

;---Load data and settings.
    version = 'v01'
    id = '2015_0416_0800'
    if n_elements(event_info) eq 0 then event_info = alfven_arc_load_data(id, event_info=event_info)

;---Plot file.
    test_plot_panels = 0
    plot_dir = event_info['plot_dir']
    plot_file = join_path([plot_dir,$
        strlowcase('fig_'+event_info.id+'_wave_'+version+'.pdf')])
    if keyword_set(test) then plot_file = 0

;---Figure out panel size.
    sgopen, 0, size=[1,1], xchsz=abs_xchsz, ychsz=abs_ychsz
    margins = [2,3,6,1]
    panel_margins = [0,0,0,0]
    uniform_ticklen = -abs_ychsz*0.15

    default_ypad = 0.4
    middle_ypad = 1
    xpads = [6,11]



end


print, fig_2015_0416_0800_wave_v01(event_info=event_info)
end