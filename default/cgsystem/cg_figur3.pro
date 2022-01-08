
pro cg_figur3
    !p.font = 1
    set_plot, 'ps'
    cgdisplay
    cgloadct, 33, rgb_table = palette
    cgplot, cgdemodata(1), layout = [2,2,3], color = 'red'
    cgcontour, cgdemodata(2), nlevel = 12, layout = [2,2,1], color = 'dodger blue'
    cgsurf, cgdemodata(2), /elevation, layout = [2,2,4], palette = palette
    cgimage, cgdemodata(19), multimargin = 4, /axes, layout = [2,2,2]
    device, /close
end