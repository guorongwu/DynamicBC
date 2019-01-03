function H = DynamicBC_connectogramOutlierA(numSlices,rextend,colorFrame2,colorFrame,mod,r2,line_widthtype,abovetype,linefactor,colortype,modcolor,nodes)
% r2, ccolor, congram.aH, th,  line_widthtype,mod,abovetype,linefactor
rextendnum = size(rextend,2);
pi_incr = 2*pi/numSlices;

Hsize = get(0,'ScreenSize');
Huse = min(Hsize(3:4))*0.8;
congram.fig = figure('pos',[Hsize(3)/2-Huse*0.5,Hsize(4)/2-Huse/2,Huse*1,Huse]);
congram.aH = axes('Parent',congram.fig,'units','norm','pos',[0.1 0.1 0.8 0.8]);
% congram.fig = figure('pos',[Hsize(3)/2-Huse*0.75,Hsize(4)/2-Huse/2,Huse*1.5,Huse]);
% congram.aH = axes('Parent',congram.fig,'units','norm','pos',[0.1*2/3 0.1 0.8*2/3 0.8]);
% congram.aH1 = axes('Parent',congram.fig,'units','norm','pos',[2/3+0.05*1/3 0.05 0.8*1/3 0.4]);
% congram.aH2 = axes('Parent',congram.fig,'units','norm','pos',[2/3+0.05*1/3 0.55 0.8*1/3 0.4]);

hold(congram.aH, 'on');
% hold(congram.aH1, 'on');
% hold(congram.aH2, 'on');
rinp = .95;
for irex = 1:rextendnum    
    if max(rextend(:,irex))==min(rextend(:,irex))
        rextendi = ones(numSlices,1);
    else
        rextendi = (rextend(:,irex)-min(rextend(:,irex)))/(max(rextend(:,irex))-min(rextend(:,irex)));
    end
    rin = 1+(irex-1)*0.2; % inner circle
    rout = 1+(irex-1)*0.2+rextendi*0.1; % outier circle
    
%     startAngle = -pi/2;
    startAngle = pi/2;
    mpX = zeros(1,numSlices);
    mpY = zeros(1,numSlices);
    th = zeros(1,numSlices);
%     midpointusedx = zeros(numSlices,1);
%     midpointusedy = zeros(numSlices,1);
    colorFrame2t = colorFrame2{irex};
    for n = 1:numSlices
        tmpC = colorFrame2t{n};
        tmpCP = colorFrame{n};
        if isempty(find(isnan(tmpCP)))
%             endAngle = startAngle+pi_incr;
%             th(n) = startAngle+pi_incr/2;
            endAngle = startAngle-pi_incr;
            th(n) = startAngle-pi_incr/2;
            t2 = linspace(startAngle, endAngle, 2);
            midang = (startAngle+endAngle)/2;
            xpoint = rinp*cos(midang);
            ypoint = rinp*sin(midang);
            startAngle = endAngle;
            
%             xin = rin*cos(t2(end:-1:1))+midpointusedx(mod(n));
%             yin = rin*sin(t2(end:-1:1))+midpointusedy(mod(n));
            xin = rin*cos(t2(end:-1:1));
            yin = rin*sin(t2(end:-1:1));
%             xout = rout(n)*cos(t2)+midpointusedx(mod(n));
%             yout = rout(n)*sin(t2)+midpointusedy(mod(n));
            xout = rout(n)*cos(t2);
            yout = rout(n)*sin(t2);
            PH(n) = patch([xout, xin],[yout, yin], tmpC, 'edgecolor', 'k', 'parent', congram.aH);
%             PH(n) = patch([xout, xin],[yout, yin], tmpC, 'edgecolor', 'k', 'parent', congram.aH1);
%             PH(n) = patch([xout, xin],[yout, yin], tmpC, 'edgecolor', 'k', 'parent', congram.aH2);
            xx = (rin + rout(n))*cos(th(n))/2;
            yy = (rin + rout(n))*sin(th(n))/2;
            %trot(n) = 180/pi* (atan((newY(2) - newY(1)) / (newX(2) - newX(1))));
            %thA=text(xx,yy,compNames{n},'horizontalalignment','center','rotation',(180/pi)*th(n),'fontsize',11,'fontweight','bold', 'color', 'k', 'tag', ['text_', num2str(n)]);
            mpX(n) = xx;
            mpY(n) = yy;
            plot(xpoint,ypoint,'parent',congram.aH,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',tmpCP,'MarkerSize',nodes);
%             plot(xpoint,ypoint,'parent',congram.aH1,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',tmpCP,'MarkerSize',nodes/2);
%             plot(xpoint,ypoint,'parent',congram.aH2,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor',tmpCP,'MarkerSize',nodes/2);
        else
%             endAngle = startAngle+pi_incr;
            endAngle = startAngle-pi_incr;
            startAngle = endAngle;
        end
    end
    
%     TH{irex} = th;
end
% line_widthtype = 1;
% abovetype = 1;
% linefactor = 1;
for i = 1:size(modcolor,1)
    ccolor(i,:) = modcolor(i,:);
end
ccolor(size(modcolor,1)+1,:) = [0.8 0.8 0.8];

DynamicBC_connectogramInner(r2, ccolor, congram.aH, th,  line_widthtype,mod,abovetype,linefactor,colortype)
% DynamicBC_schemaball2EdgeC(r2, ccolor,  congram.aH, th,  2,mod)
ROUT = ceil(max(rout)*10)/10;
axis(congram.aH,[-ROUT, ROUT, -ROUT, ROUT]*1.1)
H = congram.fig;
% DynamicBC_connectogramInnerW(r2, ccolor, congram.aH1, th,  line_widthtype,mod,abovetype,linefactor,colortype)
% DynamicBC_connectogramInnerB(r2, ccolor, congram.aH2, th,  line_widthtype,mod,abovetype,linefactor,colortype)
% axis(congram.aH1,[-ROUT, ROUT, -ROUT, ROUT]*1.1)
% axis(congram.aH2,[-ROUT, ROUT, -ROUT, ROUT]*1.1)
end
