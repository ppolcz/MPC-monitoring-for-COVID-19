function Pl_var = zf_plot_mean_var(tt,xx,xx_std,Color_Shade)
%%
%  File: plot_mean_var.m
%  Directory: 4_gyujtemegy/11_CCS/2021_COVID19_analizis/study11_SNMPC_LTV
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2021. September 10. (2021a)
%

arguments
    tt,xx,xx_std
    Color_Shade = [ 0 0.4470 0.7410 ];
end

    for sigma = 1:2

        Sh = shade(tt,xx'+sigma*xx_std',tt,xx'-sigma*xx_std','FillType',[1 2;2 1],'FillColor',Color_Shade);
        Sh(1).LineStyle = 'none';
        Sh(1).Marker = 'none';
        Sh(1).Color = Color_Shade;
        % ---
        Sh(2).LineStyle = 'none';
        Sh(2).Marker = 'none';
        Sh(2).Color = Color_Shade;
        % --- 
        Sh(3).Marker = 'none';
        Sh(3).LineStyle = 'none';
        Sh(3).FaceAlpha = 0.2/sigma;

        Pl_var(sigma) = Sh(3);

    end

    Pl_var = [ plot(tt,xx','-','Color',Color_Shade,'LineWidth',1.5) Pl_var ];

    % for alpha = 0.5:0.1:3
    % 
    %     Sh = shade(tt,xx'+alpha*xx_std',tt,xx-alpha*xx_std,'FillType',[1 2;2 1],'FillColor',Color_Shade_red);
    %     Sh(1).LineStyle = 'none';
    %     Sh(2).LineStyle = 'none';
    %     Sh(3).FaceAlpha = 0.05/alpha;
    % 
    %     if alpha == 0.5
    %         Pl_var = Sh(3);
    %     end
    % 
    % end

end
