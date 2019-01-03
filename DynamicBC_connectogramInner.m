function DynamicBC_connectogramInner(r, ccolor, gH, theta,  line_widthtype,mod,abovetype,linefactor,colortype)
sz = size(r);
t      = (0.025: 0.05 :1)';
tf        = tril(true(sz),-1);
rtf = r.*tf;
[row,col] = find(rtf);
x     = cos(theta);
y     = sin(theta);
t2  = [1-t, t].^2;
for i = 1:length(row)
    Bx = t2*[x(col(i)); x(row(i))];
    By = t2*[y(col(i)); y(row(i))];
    mod1 = mod(col(i));
    mod2 = mod(row(i));
    val = rtf(row(i),col(i));
    PlotBx{i,1} = Bx(:);
    PlotBy{i,1} = By(:);
    if mod1==mod2
        PlotLab(i,1) = 2;
        PlotCol(i,1) = mod1;
        Linevalue(i,1) = val;
        moduse(i,:) = [mod1,mod2];
    else
        PlotLab(i,1) = 1;
        PlotCol(i,1) = size(ccolor,1);
        Linevalue(i,1) = val;
        moduse(i,:) = [mod1,mod2];
    end
end

% Linevalue = (((Linevalue-min(Linevalue))/(max(Linevalue)-min(Linevalue)))*1.5)+0.5;
Linevalue = (((Linevalue-min(Linevalue))/(max(Linevalue)-min(Linevalue)))*0.5)+0.5;
Linevalue = linefactor*Linevalue;
if abovetype==1 % within above between
    [ixp iyp] = sort(PlotLab);
    for i = 1:length(row)
        if line_widthtype==1
            line_width = 2;
        else
            line_width = Linevalue(iyp(i));
        end
        if PlotLab(iyp(i),1)==2
            plot(PlotBx{iyp(i),1}, PlotBy{iyp(i),1}, 'Color',ccolor(PlotCol(iyp(i)), :) , 'LineWidth', line_width, 'parent', gH);
        else
            if colortype==1 % GRAY
                plot(PlotBx{iyp(i),1}, PlotBy{iyp(i),1}, 'Color',ccolor(PlotCol(iyp(i)), :) , 'LineWidth', line_width, 'parent', gH);
            else
                MODUSE = moduse(iyp(i),:);
                colorchanges = (ccolor(MODUSE(1),:)+ccolor(MODUSE(2),:))/2;
                plot(PlotBx{iyp(i),1}, PlotBy{iyp(i),1}, 'Color',colorchanges , 'LineWidth', line_width, 'parent', gH);
            end
        end
    end
else
    [ixp iyp] = sort(PlotLab,'descend');
    for i = 1:length(row)
        if line_widthtype==1
            line_width = 2;
        else
            line_width = Linevalue(iyp(i));
        end
        if PlotLab(iyp(i),1)==2
            plot(PlotBx{iyp(i),1}, PlotBy{iyp(i),1}, 'Color', ccolor(PlotCol(iyp(i)), :), 'LineWidth', line_width, 'parent', gH);
        else            
            if colortype==1 % GRAY
                plot(PlotBx{iyp(i),1}, PlotBy{iyp(i),1}, 'Color',ccolor(PlotCol(iyp(i)), :) , 'LineWidth', line_width, 'parent', gH);
            else
                MODUSE = moduse(iyp(i),:);
                colorchanges = (ccolor(MODUSE(1),:)+ccolor(MODUSE(2),:))/2;
                plot(PlotBx{iyp(i),1}, PlotBy{iyp(i),1}, 'Color',colorchanges , 'LineWidth', line_width, 'parent', gH);
            end
        end
    end    
end
end