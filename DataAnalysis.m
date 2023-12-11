


%% load the data
load('data');


n_f = 3:13; %indexes of the frequency bins ranging from 1Hz to 6Hz
references = {'posX', 'posY', 'excursion', 'velX', 'velY', 'vel', 'accX', 'accY', 'acc'}; %list of features
ref_analysis = references([3 6 9 4 5]);
ref_analysis{1} = 'rCOP'; ref_analysis{4} = 'vCOP_M_L'; ref_analysis{5} = 'vCOP_A_P';

conditions = {'solid_eyes_open', 'foam_eyes_open', 'solid_eyes_closed', 'foam_eyes_closed'};
% mCKC = cat(4,mCKC(:,1:2,:),mCKC(:,3:4,:));
% [p,table,stat] = CM_anovan(mCKC_BC,{'subjects' 'surf' 'ref' 'vis'});


%% topoplots 
cfg = [];
cfg.xparam = 'time';
cfg.layout = 'biosemi64.lay';
cfg.comment = 'no';
ave = [];
ref = [3:6,9]; %limited to the references 'excursion', 'vel', 'acc', 'velX' and
%'velY'


for i = 1:length(ref)
   for c = 1:length(conditions)
        ave.avg = squeeze(mean(mean(CKC_CoP(:, c, 1:63, n_f, ref(i)), 4), 1)); 
        % ave.avg = squeeze(mean(mean(CKC_CoM(:, c, 1:63, n_f, ref(i)), 4), 1)); 
        
        % cfg.zlim = [0.004 0.015];
        % cfg.zlim = [0 max(ave.avg)];
        % cfg.zlim = [0.003 0.005];
        % cfg.highlight = 'on';
        % tmp = max_elec;
        % ind = find(tmp(:, ref, 1)==0)
        % tmp(ind, :, :) = [];
        % cfg.highlightchannel = channels_min(squeeze(tmp(:, ref, 1)));
        % cfg.highlightsymbol = '*';
        % cfg.highlightcolor = [0 0 0];
        % cfg.highlightsize = 16;
        % cfg.markersize = 1;
        
        ave.fsample = 1000;
        ave.time = 0;
        ave.label = channels(1:63);
        % ave.label = Coh_COM.label(1:63);
        ave.dimord = 'chan_time';
        f=figure;  hold on;
        ft_topoplotER(cfg,ave);
        % saveas(f, strcat('CKCFigures/Traces2/COP_', conditions{c}, '_', references{i}, '_V4.svg'));
        %         close(f);
   end
end

%% number of participants showing significant CKC
%dimensions for sig and p: participants*conditions*references
sig_cop = p_CKC_CoP < 0.05/(length(n_f)*length(SM1));
sig_com = p_CKC_CoM < 0.05/(length(n_f)*length(SM1));
%dimensions for num and percentage: conditions*references
num_part_cop = squeeze(sum(sig_cop, 1));
num_part_com = squeeze(sum(sig_cop, 1));
percentage_cop = num_part_cop./36*100;
percentage_com = num_part_com./36*100;


%% Effect of the condition on the number of available epochs for CKC computation
[p, tbl, stat] = anovan(reshape(NAVE, 36*4, 1), {indiv, cond}, 'varnames', {'participant', 'condition'});

%% violin plots for CKC COP & COM
% separate plots for COP and COM
fig = figure('Units','centimeters','Position',[0 0 16 20])

ref_order = [6 5 3 4 9];
for n_center = 1:2
    if n_center == 1
        this_CKC = maxCKC_CoP(:, :, ref_order);
        center_type = 'COP';
    else
        this_CKC = maxCKC_CoM(:, :, ref_order);
        center_type = 'COM';
    end
    for n_ref = 1:5
        ymax = 0.25;
        if n_ref >3 | n_center == 2
            ymax = 0.15;
        end
        % subplot(1,5,n_ref)
        xx = 1+(mod(n_ref-1,3))*5;
        yy = 5 - 3*(n_ref > 3);
        if n_center == 1;
            yy = yy + 7;
        end
        axes('Units','centimeters','Position',[xx yy 5 15*ymax],'FontSize',7,'FontName','Arial');
        CM_violin_plot_TL(atanh(sqrt(this_CKC(:,:,n_ref))),0)
        title(references{ref_order(n_ref)})
        ylim([0 ymax])
        xlim([0 15])
        xticks([])
        if any(n_ref == [2 3 5])
            yticks([])
        end
    end
end
% print('-dsvg',fullfile('/home/mathieu/Documents/Articles/publi/Thomas_CKC_COP/data/',['CKC_CoP_CoM_violin.svg']))
% close(fig)


% grouped plot for COP and COM
fig = figure('Units','centimeters','Position',[0 0 16 12])
ref_order = [2 5 1 4 3];
this_CKC = cat(1, maxCKC_CoP(:, :, ref), maxCKC_CoM(:, :, ref));
for n_ref = 1:5
    ymax = 0.25;
    if n_ref > 3
        ymax = 0.15;
    end
    xx = 1+(mod(n_ref-1,3))*5;
    yy = 5 - 4*(n_ref > 3);
    axes('Units','centimeters','Position',[xx yy 5 20*ymax],'FontSize',7,'FontName','Arial');
    CM_violin_plot_TL(atanh(sqrt(this_CKC(:,:,ref_order(n_ref)))),0,size(this_CKC,1)/2)
    title(references{ref_order(n_ref)})
    ylim([0 ymax])
    xlim([0 15])
    xticks([])
    if any(n_ref == [2 3 5])
        yticks([])
    end
end
% print('-dsvg',fullfile('/home/mathieu/Documents/Articles/publi/Thomas_CKC_COP/data/CKC_CoP_CoM_violin_grouped.svg'))
% close(fig)

%% ANOVA for CoP and CoM average speed, antero-posterior and medio-lateral standard deviation

%%%CoP%%%
[p,table,stat] = CM_anovan(average_speed_CoP,{'subjects' 'cond'});
table(:,[1:3 11 6 7])
[h,p,tmp,stat] = ttest(average_speed_CoP(:,[1 1 1 2 2 3])-average_speed_CoP(:,[2 3 4 3 4 4]))

[p,table,stat] = CM_anovan(std_CoPX,{'subjects' 'cond'});
table(:,[1:3 11 6 7])
[h,p,tmp,stat] = ttest(std_CoPX(:,[1 1 1 2 2 3])-std_CoPX(:,[2 3 4 3 4 4]))

[p,table,stat] = CM_anovan(std_CoPY,{'subjects' 'cond'});
table(:,[1:3 11 6 7])
[h,p,tmp,stat] = ttest(std_CoPY(:,[1 1 1 2 2 3])-std_CoPY(:,[2 3 4 3 4 4]))

figure; b = violin(average_speed_CoP*100*1000, 'xlabel', conditions, 'facecolor', [0,1,1; 0,0,1; 0.8,0.1,0.9; 0.5,0.2,0.6]);
figure; b = violin(std_CoPX*100, 'xlabel', conditions, 'facecolor', [0,1,1; 0,0,1; 0.8,0.1,0.9; 0.5,0.2,0.6]);
figure; b = violin(std_CoPY*100, 'xlabel', conditions, 'facecolor', [0,1,1; 0,0,1; 0.8,0.1,0.9; 0.5,0.2,0.6]);

%%%CoM%%%
[p,table,stat] = CM_anovan(average_speed_CoM,{'subjects' 'cond'});
table(:,[1:3 11 6 7])
[h,p,tmp,stat] = ttest(average_speed_CoM(:,[1 1 1 2 2 3])-average_speed_CoM(:,[2 3 4 3 4 4]))

[p,table,stat] = CM_anovan(std_CoMX,{'subjects' 'cond'});
table(:,[1:3 11 6 7])
[h,p,tmp,stat] = ttest(std_CoMX(:,[1 1 1 2 2 3])-std_CoMX(:,[2 3 4 3 4 4]))

[p,table,stat] = CM_anovan(std_CoMY,{'subjects' 'cond'});
table(:,[1:3 11 6 7])
[h,p,tmp,stat] = ttest(std_CoMY(:,[1 1 1 2 2 3])-std_CoMY(:,[2 3 4 3 4 4]))

figure; b = violin(average_speed_CoM*100*1000, 'xlabel', conditions, 'facecolor', [0,1,1; 0,0,1; 0.8,0.1,0.9; 0.5,0.2,0.6]);
figure; b = violin(std_CoMX*100, 'xlabel', conditions, 'facecolor', [0,1,1; 0,0,1; 0.8,0.1,0.9; 0.5,0.2,0.6]);
figure; b = violin(std_CoMY*100, 'xlabel', conditions, 'facecolor', [0,1,1; 0,0,1; 0.8,0.1,0.9; 0.5,0.2,0.6]);

%% ANOVA for CKC one reference at a time
for n_center = 1:2 % 1: CoP; 2: CoM
    for n_ref = 1:9
        if n_center == 1
            mCKC_BC = maxCKC_CoP(:,:,n_ref);
        else
            mCKC_BC = maxCKC_CoM(:,:,n_ref);
        end
        mCKC_BC(:) = boxcox(mCKC_BC(:));
        [p,table,stat] = CM_anovan(mCKC_BC,{'subjects' 'cond'});
        references{n_ref}
        table(:,[1:3 11 6 7])
        [h,p,tmp,stat] = ttest(mCKC_BC(:,[1 1 1 2 2 3])-mCKC_BC(:,[2 3 4 3 4 4]))
    end
end

% significant difference between aCOP and vCOPml
[h,p,tmp,stat] = ttest(maxCKC_CoP(:,:,4)-maxCKC_CoP(:,:,9))


n_ref = 3;
Nsim = 1000;
Nsub = size(maxCKC_CoP,1);
p_dist = zeros(2,Nsim);
for n_sim = 1:Nsim
    keep_sub = ceil(rand(1,Nsub)*Nsub);
    for n_center = 1:2
        if n_center == 1
            mCKC_BC = maxCKC_CoP(keep_sub,:,n_ref);
        else
            mCKC_BC = maxCKC_CoM(keep_sub,:,n_ref);
        end
        mCKC_BC(:) = boxcox(mCKC_BC(:));
        [p,table,stat] = CM_anovan(mCKC_BC,{'subjects' 'cond'});
        p_dist(n_center,n_sim) = table{3,[7]};
    end
end
p = [mean(diff(p_dist) > 0) mean(diff(p_dist) < 0)];
p = max(min(p,1-p))*2


% Difference between CKC for COP vs. COM
[h,p,tmp,stat] = ttest(maxCKC_CoP(:, :, ref)-maxCKC_CoM(:, :, ref))
references
squeeze(stat.tstat)

%% Comparison between spectral profiles using a permutation test

% % % % Comparison with permutation statistics. Not optimal for kurtosis 
% 
n_ref = [5 6];
Nsim = 1000;
Nsub = size(CKC_FEC,1);
p_dist = zeros(2,Nsim);

spctr = squeeze(mean(CKC_FEC(:,2:21,n_ref),1));
% spctr = spctr./sqrt(sum(spctr.^2));
spctr = spctr-mean(spctr(16:20,:));
% d = sum(diff(spctr,[],2).^2);
d = diff(kurtosis(spctr));

d_dist = zeros(1,Nsim);
for n_sim = 1:Nsim
    flip = find(rand(1,Nsub)>0.5);
    not_flip = 1:Nsub;
    not_flip(flip) = [];
    spctr = squeeze(mean(cat(1,CKC_FEC(not_flip,2:21,n_ref),CKC_FEC(flip,2:21,n_ref([2 1]))),1));
    % spctr = spctr./sqrt(sum(spctr.^2));
    spctr = spctr-mean(spctr(16:20,:));
    % d_dist(n_sim) = sum(diff(spctr,[],2).^2);
    d_dist(n_sim) = diff(kurtosis(spctr));
end
p = mean(d <= d_dist)
% 


n_ref = [5 6];
Nsub = size(CKC_FEC,1);

spctr = squeeze(mean(CKC_FEC(:,2:21,n_ref),1));
% param = kurtosis(spctr)
% param = median(spctr);
% param = find(cumsum(spctr) <= sum(spctr)/2, 2, 'last')

Nsim = 10000;
param_dist = zeros(2,Nsim);
for n_sim = 1:Nsim
    keep_sub = ceil(rand(1,Nsub)*Nsub);
    spctr = squeeze(mean(CKC_FEC(keep_sub,2:21,n_ref),1));
    param_dist(:,n_sim) = kurtosis(spctr);
%     param_dist(1,n_sim) = find(cumsum(spctr(:,1)) <= sum(spctr(:,1))/2, 1, 'last');
%     param_dist(2,n_sim) = find(cumsum(spctr(:,2)) <= sum(spctr(:,2))/2, 1, 'last');
end
p = [mean(diff(param_dist) > 0) mean(diff(param_dist) < 0)];
p = max(min(p,1-p))*2


CKC_FEC_sig_5 = CKC_FEC(sig(:,4,5),2:21,5);
CKC_FEC_sig_6 = CKC_FEC(sig(:,4,6),2:21,6);
Nsub5 = size(CKC_FEC_sig_5, 1); Nsub6 = size(CKC_FEC_sig_6, 1);
Nsim = 10000;
param_dist = zeros(2,Nsim);
for n_sim = 1:Nsim
    keep_sub5 = ceil(rand(1,Nsub5)*Nsub5);
    keep_sub6 = ceil(rand(1,Nsub6)*Nsub6);
    spctr5 = squeeze(mean(CKC_FEC_sig_5(keep_sub5,:,1),1));
    spctr6 = squeeze(mean(CKC_FEC_sig_6(keep_sub6,:,1),1));

    param_dist(1,n_sim) = kurtosis(spctr5);
    param_dist(2,n_sim) = kurtosis(spctr6);
%     param_dist(1,n_sim) = find(cumsum(spctr(:,1)) <= sum(spctr(:,1))/2, 1, 'last');
%     param_dist(2,n_sim) = find(cumsum(spctr(:,2)) <= sum(spctr(:,2))/2, 1, 'last');
end
p = [mean(diff(param_dist) > 0) mean(diff(param_dist) < 0)];
p = max(min(p,1-p))*2


%% Delays


figure; violin({delays_yfm_rCOP, delays_yfm_vCOPAP, delays_yfm_vCOP}, 'xlabel', {'rCOP', 'vCOP_{AP}', 'vCOP'}, 'facecolor', [0.5,0.2,0.6;0.5,0.2,0.6;0.5,0.2,0.6]);
figure; hold on;
fields = fieldnames(phases.rCOP);
for i = 1:length(fields)
    arr = phases.rCOP.(fields{i});
    mdl = fitlm(arr(:,1), arr(:,2));
    plot(arr(:,1), arr(:,2), 'k--o', 'LineWidth', 1);
    plot(arr(:,1), arr(:,1)*mdl.Coefficients.Estimate(2)+mdl.Coefficients.Estimate(1), '-','Color', [0.5, 0.2, 0.6], 'LineWidth', 2);
end
xlim([0,10]); ylim([-180,180]); yticks([-180, -90, 0, 90, 180]);
xlabel('Frequency (Hz)'); ylabel('Angle (°)');

figure; hold on;
fields = fieldnames(phases.vCOPAP);
for i = 1:length(fields)
    arr = phases.vCOPAP.(fields{i});
    mdl = fitlm(arr(:,1), arr(:,2));
    plot(arr(:,1), arr(:,2), 'k--o', 'LineWidth', 1);
    plot(arr(:,1), arr(:,1)*mdl.Coefficients.Estimate(2)+mdl.Coefficients.Estimate(1), '-','Color', [0.5, 0.2, 0.6], 'LineWidth', 2);
end
xlim([0,10]); ylim([-180,180]); yticks([-180, -90, 0, 90, 180]);
xlabel('Frequency (Hz)'); ylabel('Angle (°)');

figure; hold on;
fields = fieldnames(phases.vCOP);
for i = 1:length(fields)
    arr = phases.vCOP.(fields{i});
    mdl = fitlm(arr(:,1), arr(:,2));
    plot(arr(:,1), arr(:,2), 'k--o', 'LineWidth', 1);
    plot(arr(:,1), arr(:,1)*mdl.Coefficients.Estimate(2)+mdl.Coefficients.Estimate(1), '-','Color', [0.5, 0.2, 0.6], 'LineWidth', 2);
end
xlim([0,10]); ylim([-180,180]); yticks([-180, -90, 0, 90, 180]);
xlabel('Frequency (Hz)'); ylabel('Angle (°)');

%Wilcoxon signed rank test to determine if the estimated delays are
%significantly different from 0
signrank(delays_yfm_rCOP)
signrank(delays_yfm_vCOP)
signrank(delays_yfm_vCOPAP)


%% Behavioral relevance

n_ref = [3 6 5];
x1 = (std_CoPY(:,4)-std_CoPY(:,1))./(std_CoPY(:,4)+std_CoPY(:,1));
x2 = (average_speed_CoP(:,4)-average_speed_CoP(:,1))./(average_speed_CoP(:,4)+average_speed_CoP(:,1));
x3 = std_CoPY(:,4)*1000; % in mm
x4 = average_speed_CoP(:,4)*1e5; % in cm/s
y1 = squeeze((maxCKC_CoP(:,4,n_ref)-maxCKC_CoP(:,1,n_ref))./(maxCKC_CoP(:,4,n_ref)+maxCKC_CoP(:,1,n_ref)));
% y2 = squeeze((mCKC(:,4,n_ref)+mCKC(:,1,n_ref)-mCKC(:,2,n_ref)-mCKC(:,3,n_ref))./sum(mCKC(:,:,n_ref),2));
y2 = squeeze(atanh(sqrt(maxCKC_CoP(:,4,n_ref))));

msize = 4;
Nsim = 1000;


%CCA between relative increase in CoP instability parameters (average speed, AP SD) 
% and relative increase in CKC (excurion, vCoP_AP, vCoP) between quiet
% standing and standing on foam eyes closed
x = cat(2,x1,x2);
y = y1;

stat = CM_CCA_leave_one_out(x',y',Nsim)

fig = figure('Units','centimeters','Position',[0 0 25 20])
axes('Units','centimeters','Position',[1 1 4 3],'FontSize',7,'FontName','Arial');
for n_x = 1:length(stat.Vx)
    plot(n_x+[-0.4 0.4 0.4 -0.4 -0.4],stat.Vx(n_x)*[0 0 1 1 0],'k');
    hold on
end
for n_y = 1:length(stat.Vy)
    plot(length(stat.Vx)+1+n_y+[-0.4 0.4 0.4 -0.4 -0.4],stat.Vy(n_y)*[0 0 1 1 0],'k');
end
plot([0 length(stat.Vx)+length(stat.Vy)+2],[0 0],'r')
ylim([-1 1])

axes('Units','centimeters','Position',[7 1 4 3],'FontSize',7,'FontName','Arial');
xx = stat.x; yy = stat.y;
plot(xx,yy,'ok','markersize',msize)
pp = polyfit(xx,yy,1);
xx = [min(xx),max(xx)]; yy = pp(2)+pp(1)*xx;
hold on
plot(xx,yy,'r')
axis([-2.1 2.1 -4 2.5])
title(['r = ' num2str(stat.r,3) ', p = ' num2str(stat.p,3)])
xlabel('CoP canonical variate')
ylabel('CKC canonical variate')
set(gca,'fontsize',7)

% stat = CM_pls(x',y',100)

% xx = (x-mean(x,1))./std(x,[],1);
% yy = (y-mean(y,1))./std(y,[],1);
% yy = y;
axes('Units','centimeters','Position',[13 1 4 3],'FontSize',7,'FontName','Arial');
xx = mean(x,2); yy = mean(y,2);
plot(xx,yy,'ok','markersize',msize)
[r,p] = corr(xx,yy);
pp = polyfit(xx,yy,1);
xx = [min(xx),max(xx)]; yy = pp(2)+pp(1)*xx;
hold on
plot(xx,yy,'r')
title(['r = ' num2str(r,3) ', p = ' num2str(p,3)])
axis([0.15 0.45 -0.2 1])
xlabel('Increase in instability')
ylabel('Increase in CKC')
set(gca,'fontsize',7)

% CCA between absolute CoP instability parameters (average speed, AP SD) 
% and absolute CKC (excurion, vCoP_AP, vCoP) when standing on foam eyes
% closed
x = cat(2,x3,x4);
y = y2;

stat = CM_CCA_leave_one_out(x',y',Nsim)

axes('Units','centimeters','Position',[1 6 4 3],'FontSize',7,'FontName','Arial');
for n_x = 1:length(stat.Vx)
    plot(n_x+[-0.4 0.4 0.4 -0.4 -0.4],stat.Vx(n_x)*[0 0 1 1 0],'k');
    hold on
end
for n_y = 1:length(stat.Vy)
    plot(length(stat.Vx)+1+n_y+[-0.4 0.4 0.4 -0.4 -0.4],stat.Vy(n_y)*[0 0 1 1 0],'k');
end
plot([0 length(stat.Vx)+length(stat.Vy)+2],[0 0],'r')
ylim([-1 1])

axes('Units','centimeters','Position',[7 6 4 3],'FontSize',7,'FontName','Arial');
xx = stat.x; yy = stat.y;
plot(xx,yy,'ok','markersize',msize)
pp = polyfit(xx,yy,1);
xx = [min(xx),max(xx)]; yy = pp(2)+pp(1)*xx;
hold on
plot(xx,yy,'r')
axis([-2.3 2.6 -1 1.5])
title(['r = ' num2str(stat.r,3) ', p = ' num2str(stat.p,3)])
xlabel('CoP canonical variate')
ylabel('CKC canonical variate')
set(gca,'fontsize',7)


axes('Units','centimeters','Position',[13 6 4 3],'FontSize',7,'FontName','Arial');
% % % x3 = std_COPY(:,4)*1000; % en mm
% % % x4 = COP_vel(:,4)*1e5; % en cm/s
xx = x(:,2)-x(:,1)/3; yy = y(:,2)/2-y(:,1); % xx: vCOP`/ 3 cm/s     stdCOP / 3 mm
plot(xx,yy,'ok','markersize',msize)
[r,p] = corr(xx,yy);
pp = polyfit(xx,yy,1);
xx = [min(xx),max(xx)]; yy = pp(2)+pp(1)*xx;
hold on
plot(xx,yy,'r')
title(['r = ' num2str(r,3) ', p = ' num2str(p,3)])
axis([0 9 -0.05 0.03])
xlabel('Speed(COP_{AP})/(3 cm/s) - SD(COP_{AP})/(1 mm)')
ylabel('CKC(vCOP)/2 - CKC(rCOP)')
set(gca,'fontsize',7)


axes('Units','centimeters','Position',[19 6 4 3],'FontSize',7,'FontName','Arial');
xx = x(:,2)+x(:,1)/3; yy = mean(y,2);
plot(xx,yy,'ok','markersize',msize)
[r,p] = corr(xx,yy);
pp = polyfit(xx,yy,1);
xx = [min(xx),max(xx)]; yy = pp(2)+pp(1)*xx;
hold on
plot(xx,yy,'r')
title(['r = ' num2str(r,3) ', p = ' num2str(p,3)])
axis([6 15 0.04 0.18])
xlabel('Speed(COP_{AP})/(3 cm/s) + SD(COP_{AP})/(1 mm)')
ylabel('global CKC for CoP')
set(gca,'fontsize',7)



%CCA between relative increase in CoP instability parameters (average
%speed, AP SD) between quiet stand and standing on foam eyes closed
% and absolute CKC (excurion, vCoP_AP, vCoP) when standing on foam eyes closed
x = cat(2,x1,x2);
y = y2;

stat = CM_CCA_leave_one_out(x',y',Nsim)

axes('Units','centimeters','Position',[1 11 4 3],'FontSize',7,'FontName','Arial');
for n_x = 1:length(stat.Vx)
    plot(n_x+[-0.4 0.4 0.4 -0.4 -0.4],stat.Vx(n_x)*[0 0 1 1 0],'k');
    hold on
end
for n_y = 1:length(stat.Vy)
    plot(length(stat.Vx)+1+n_y+[-0.4 0.4 0.4 -0.4 -0.4],stat.Vy(n_y)*[0 0 1 1 0],'k');
end
plot([0 length(stat.Vx)+length(stat.Vy)+2],[0 0],'r')
ylim([-1 1])

axes('Units','centimeters','Position',[7 11 4 3],'FontSize',7,'FontName','Arial');
xx = stat.x; yy = stat.y;
plot(xx,yy,'ok','markersize',msize)
pp = polyfit(xx,yy,1);
xx = [min(xx),max(xx)]; yy = pp(2)+pp(1)*xx;
hold on
plot(xx,yy,'r')
axis([-2 2.3 -2 5])
title(['r = ' num2str(stat.r,3) ', p = ' num2str(stat.p,3)])
xlabel('CoP canonical variate')
ylabel('CKC canonical variate')
set(gca,'fontsize',7)

% stat = CM_pls(x',y',100)

% xx = (x-mean(x,1))./std(x,[],1);
% yy = (y-mean(y,1))./std(y,[],1);
% yy = y;
axes('Units','centimeters','Position',[13 11 4 3],'FontSize',7,'FontName','Arial');
xx = mean(x,2); yy = mean(y,2);
plot(xx,yy,'ok','markersize',msize)
[r,p] = corr(xx,yy);
pp = polyfit(xx,yy,1);
xx = [min(xx),max(xx)]; yy = pp(2)+pp(1)*xx;
hold on
plot(xx,yy,'r')
title(['r = ' num2str(r,3) ', p = ' num2str(p,3)])
axis([0.15 0.45 0.04 0.18])
xlabel('Increase in instability')
ylabel('global CKC for CoP')
set(gca,'fontsize',7)




    