function fig_pos=get_fig_pos(fig_num,fig_rows)

%fig
set(0,'Units','pixels') 
screen_size=get(0,'ScreenSize');
screen_height=screen_size(4);
screen_width=screen_size(3);
set(0,'Units','inches') 
screen_size_i=get(0,'ScreenSize');
res=screen_size(4)/screen_size_i(4);

set(0,'Units','pixels')
fig_per_row=ceil(fig_num/fig_rows);

fig_height=(screen_height-1.25*res)/fig_rows;


fig_left=0;
for k=1:fig_num
    fig_width=screen_width/fig_per_row;
    if k <= fig_per_row
        fig_botom=screen_height-fig_height;
        fig_pos{k}=[floor(fig_left),floor(fig_botom),floor(fig_width),floor(fig_height)];
        fig_left=fig_left+fig_width;
    else
        if k== fig_per_row+1
        fig_left=0;
        end
        fig_botom=screen_height-(ceil(k/fig_per_row))*(fig_height-.25*res);
        fig_pos{k}=[floor(fig_left),floor(fig_botom),floor(fig_width),floor(fig_height)];
        fig_left=fig_left+fig_width;
    end
    
    
end

set(0,'DefaultAxesFontName','Times')
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultTextFontName','Times')
set(0,'DefaultTextFontSize',18)
set(0,'DefaultLineLineWidth',1.1)
