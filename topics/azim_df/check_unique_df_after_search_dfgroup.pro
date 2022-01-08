;+
; Check how many unique vertex are found after divided into subgroups.
;-


    events = azim_df_find_dfgroup(project=project)
    ;events = azim_df_find_subgroup(project=project)
    lprmsg, '# of events: '+string(events.length,format='(I0)')
    ndf = 0
    foreach event, events do ndf += event.df_list.length
    lprmsg, '# of DFs: '+string(ndf,format='(I0)')
end