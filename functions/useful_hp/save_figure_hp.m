function save_figure_hp(h,save_flag_png,png_outfile,dpi,save_flag_fig,fig_outfile)

if save_flag_png
    set(h,'PaperPositionMode','auto')
    print(h,png_outfile,'-dpng',['-r',dpi])
end

if save_flag_fig
    savefig(h,fig_outfile)
end

end