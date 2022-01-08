;+
; Read all substorm_list.
;-

root_dir = join_path([googledir(),'works','works','azim_dp_vs_fac'])
file = join_path([root_dir,'data','azim_dp_vs_fac_substorm_list.txt'])
plot_dir = join_path([root_dir,'plot'])
candidate_list = azim_dp_candidate_read(filename=file)
foreach candidate, candidate_list do begin
    time_range = candidate.time_range+[-1,1]*1800d
    probes = candidate.probes
    mlt_range = candidate.mlt_range
    plot_file = join_path([plot_dir,'fig_ewogram_of_dp_and_up_down_current_'+strjoin(time_string(time_range,tformat='YYYY_MMDD_hh'),'_')+'_v01.pdf'])
    fig_ewogram_of_dp_and_up_down_current, time_range, probes=probes, mlt_range=mlt_range, filename=plot_file
endforeach

end
