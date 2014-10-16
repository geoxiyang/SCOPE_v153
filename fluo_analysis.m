%% fluo_analysis

% Author: Xi Yang (geoxiyang@gmail.com)

% Date: Aug. 8, 2014

% Read the SCOPE-produced fluorescence file
% Plot fluorescence with one parameter changes and the rest fixed
% Calculate the ratio between fluorescence at a single wavelength and 
% the integral of SIF

% In the future, this file make a LUT, and use it as the input to SCOPE,
% and then we plot the response of fluorescence as a function of variables

clear all
clc

global constants
[constants] = define_constants();


work_dir    = '/Volumes/XiYangResearch/src/SCOPE/SCOPE_v1.53/output/HF_ts_2014-09-17-1338/';      % Change the input folder for the run you want to plot
ofigure     = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/';
time_dir    = '/Volumes/XiYangResearch/src/SCOPE/data/input/dataset HF_ts/';

wl          = dlmread([work_dir 'wl.dat'],'',2,0);
fEnergy     = dlmread([work_dir 'fluorescence.dat'],'',2,0);
time        = dlmread([time_dir 't_.dat'],'',0,0);

% f740        = fEnergy(:,740-639);
% f737        = fEnergy(:,737-639);
 f760        = fEnergy(:,760-639);                     % Common O2A retrieval wavelength
% f755        = fEnergy(:,755-639);                     % at the GOSAT wavelenght. Unit is: W m-2 s-1 um-1 sr-1

% plot(time(1:4120,1),f760,'bo')
sif_day = zeros(254-170+1,1);

for ii = 1:254-170
    lb = ii-1+170;
    ub = ii+170;
    
    time_sub = time>=lb & time <ub;
    sif_day(ii,1)  = nanmean(f760(time_sub));
    
    
    
    
end
plot(170:254,sif_day,'bo')





% for ii = 1:size(f760)
% %     fPhotons    = 1E6   * e2phot(1E-9*wl(wl>639 & wl<851),pi*fEnergy(ii,:));   % convert into umol m-2 s-1 um-1
% %     fPhotonsInt = 0.001 * Sint(fPhotons,wl(wl>639 & wl<851));    % convert into umol m-2 s-1
% 
%     Factor(ii)      = f755(ii)/sum(fEnergy(ii,:));
% %     Factor2(ii)     = f760(ii)/f737(ii);
% %     Factor3(ii)     = f755(ii)/f737(ii);
% %     Factor4(ii)     = f760(ii)/f740(ii);
% end


% para_combination = dlmread([work_dir 'pars_and_input_short.dat'],'',1,0);
% plot(para_combination,Factor','r*-');
% 
% xlabel('V_{cmax}(umol m^{-2} s^{-1})','FontSize',16,'FontName','Helvetica');
% ylabel('F755/SIF','FontSize',16,'FontName','Helvetica');
% set(gca,'FontSize',16,'FontName','Helvetica');
% print('-dpng','-r300',[ofigure 'factor-vcmax-chl40.png']);
% print('-depsc','-r300',[ofigure 'factor-vcmax-chl40.eps']);


%% Linear regression
% Between Factors and independent varibles vcmax etc.
% x = [ones(size(para_combination(:,1))),para_combination];
% 
% [b,bint,r,rint,stats] = regress(Factor',x);
% 
% b
% stats
% 
% plot(b'*x',Factor,'bo',min(Factor):(max(Factor)-min(Factor))/100.0:max(Factor),min(Factor):(max(Factor)-min(Factor))/100.0:max(Factor),'r-')
% xlabel('Regression F755/SIF','FontSize',16,'FontName','Helvetica');
% ylabel('SCOPE F755/SIF','FontSize',16,'FontName','Helvetica');
% set(gca,'FontSize',16,'FontName','Helvetica');
% print('-dpng','-r300',[ofigure 'factor-vcmax-chl-regress_newSCOPE.png']);
% print('-depsc','-r300',[ofigure 'factor-vcmax-chl-regress_newSCOPE.eps']);



%% Draw figures
% rawdata = [Factor',para_combination];
% [pi,pj] = size(rawdata);
% 
% uni_para = unique(rawdata(:,2));
% fig = figure;
% cc  = hsv(10);
% M = zeros(size(uni_para));
% 
% xlabel('V_{cmax}(umol m^{-2} s^{-1})','FontSize',16,'FontName','Helvetica');
% ylabel('F755/SIF','FontSize',16,'FontName','Helvetica');
% xlim([0,130]);
% set(gca,'FontSize',16,'FontName','Helvetica');
% 
% for jj = 1:size(uni_para)
%     ind  = find(rawdata(:,2) == uni_para(jj));
%     data2 = rawdata(ind,:);
%     hold on
%     plot(data2(:,3),data2(:,1),'-*','Color',cc(jj,:));
%     M(jj) = uni_para(jj);
% end
% hold off
% outM = strtrim(cellstr(num2str(M))');
% 
% legend(gca,outM{:});
% print(gcf,'-dpng','-r300',[ofigure 'factor-vcmax-chl_newSCOPE.png']);
% print(gcf,'-depsc','-r300',[ofigure 'factor-vcmax-chl_newSCOPE.eps']);                                          
