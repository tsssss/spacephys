;+
; Check detected DFs with long width.
;-


search_step = 'df'
dfs = list()
routines = 'azim_df_search_'+['pre','post']+'_midn_events'
foreach routine, routines do begin
    events = call_function(routine, search_step=search_step)
    foreach event, events do dfs.add, event.df_list, /extract
endforeach
ndf = dfs.length
widths = fltarr(ndf)
foreach df, dfs, df_id do widths[df_id] = df.width
heights = fltarr(ndf)
foreach df, dfs, df_id do heights[df_id] = df.height


ntarget_eg = 50
project = azim_df_load_project()
data_type = 'theta'

min_width = 900.
index = where(widths ge min_width, count)

min_height = 60
index = where(heights ge min_height, count)


if count ne 0 then begin
    the_dfs = dfs[index]
    if count ge ntarget_eg then begin
        index = floor(findgen(ntarget_eg)*(count-1)/(ntarget_eg-1))
        the_dfs = the_dfs[index]
    endif
    
    foreach df, the_dfs do begin
        probe = df.probe
        df_time_range = df.time_range
        time_range = df_time_range+[-1,1]*width
        azim_df_read_data, data_type, probe=probe, time_range=time_range, project=project

        sgopen, 0, xsize=4, ysize=3
        tpos = sgcalcpos(1, xchsz=xchsz, ychsz=ychsz)
        the_var = probe+'_'+data_type
        options, the_var, 'constant', 0
        tplot, the_var, trange=time_range, position=tpos
        timebar, df_time_range, linestyle=1
        tx = tpos[0]
        ty = tpos[3]+ychsz*0.5
        width = df.width
        height = df.height
        xyouts, tx,ty,/normal, strupcase(probe)+' DF, width (min): '+$
            string(width/60,format='(I0)')+', height (deg): '+$
            string(height,format='(F5.1)')


        stop
    endforeach
endif

end