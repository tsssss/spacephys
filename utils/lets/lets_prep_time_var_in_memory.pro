function lets_prep_time_var_in_memory, time_range, time_step, output=time_var

    times = make_bins(time_range, time_step)
    store_data, time_var, times, times, limits={time_step:time_step, time_range:time_range}
    return, time_var

end