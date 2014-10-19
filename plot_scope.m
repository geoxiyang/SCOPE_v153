%% Load fluorescence data and plot them

clear all
clc

%% 1. Load data

work_dir    = '/Volumes/XiYangResearch/src/SCOPE/SCOPE_v1.53/output/HF_2013_2014-10-16-1036/';
wl          = dlmread([work_dir 'wl.dat'],'',2,0);
fEnergy     = dlmread([work_dir 'fluorescence.dat'],'',2,0);
A_fluxes    = dlmread([work_dir 'fluxes.dat'],'',[2,10,49,10]);
Rin         = dlmread([work_dir 'pars_and_input_short.dat'],'',1,0);

wl          = wl(1,wl(1,:) >=640 & wl(1,:)<=850);
scope_f760  = fEnergy(:,wl(1,:) == 760);
PAR         = Rin * 4.5;

% plot(PAR,scope_f760,'r-')
% plot(PAR,A_fluxes,'b-')
% plot(scope_f760,A_fluxes,'k-','LineWidth',1.5);
% [A,H1,H2] = plotyy(PAR,scope_f760,PAR,A_fluxes);

% ofigure     = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/GPP_SIF_PAR.png';
% xlabel('SIF(mw m^{-2} sec^{-1} sr^{-1})','FontSize',20,'FontName','Whitney-Book');
% ylabel('GPP(umol m^{-2} sec^{-1})','FontSize',20,'FontName','Whitney-Book');
% set(H1,'color','red');
% set(H2,'color','blue');
% set(A(1),'YColor','r','FontSize',35,'FontName','Whitney-Book');
% set(A(2),'YColor','b','FontSize',35,'FontName','Whitney-Book');
% xlabel('PPFD(umol m^{-2} sec^{-1})');
% ylabel(A(1),'SIF(mw m^{-2} sec^{-1} sr^{-1})');
% ylabel(A(2),'GPP(umol m^{-2} sec^{-1})');
% a = get(gcf,'OuterPosition');

% set(gcf,'Units','normalized','Position',[0.2,0.2,0.4,1],'OuterPosition',[0,0,1,1]);
%set(gcf,'Units','normalized','OuterPosition',[0,0,1,1]);
%set(gcf,'OuterPosition',[a(1),a(2),a(3),a(4)]);

% set(gca,'FontSize',35,'FontName','Whitney-Book','Position',[0.13 0.2 0.5 0.7])
% set(gcf,'paperPositionMode','auto') % make the print as big as the figure
% print(gcf, '-dpng','-r300', ofigure);

%legend('Modelled GPP','Measured GPP');

load('/Volumes/XiYangResearch/src/HF_Fluo_data/SIF760daily.mat','halfhourly_result');
load('/Volumes/XiYangResearch/src/HF_Fluo_data/HF_2013_GPP.mat','gpp_raw')


%% 2. Plot data

ofigure     = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/SIF_test6.png';
ofigure1    = '/Volumes/XiYangResearch/Projects/1.SCOPE_HF/1.JPG/gpp_test6.png';


% 2.1 SIF obs vs. model
wl          = wl(1,wl(1,:) >=640 & wl(1,:)<=850);
scope_f760  = fEnergy(:,wl(1,:) == 760);

day         = 175; 
measure_f760= halfhourly_result(floor(halfhourly_result(:,1)) == day,:);

scope_f760 = [measure_f760(:,1),scope_f760];

figure
plot(scope_f760(:,1)-floor(scope_f760(:,1)),scope_f760(:,2),'bo',measure_f760(:,1)-floor(measure_f760(:,1)),measure_f760(:,2),'r*')
xlabel('Hour','FontSize',20,'FontName','Helvetica');
ylabel('SIF(w m-2 nm-1 sr-1)','FontSize',20,'FontName','Helvetica');
set(gca,'FontSize',20,'FontName','Helvetica')
legend('Modelled SIF','Measured SIF');

set(gcf,'paperPositionMode','auto') % make the print as big as the figure
print(gcf, '-dpng','-r300', ofigure);

% 2.2 GPP obs vs. model
day         = 175;
hours       = (halfhourly_result(floor(halfhourly_result(:,1)) == day,1)-day)*24.0;
measure_gpp = gpp_raw(floor(halfhourly_result(:,1)) == day,:);

figure
plot(hours,A_fluxes,'bo',hours,measure_gpp,'r*')

xlabel('Hour','FontSize',20,'FontName','Whitney-Book');
ylabel('GPP(umol m^{-2} s^{-1})','FontSize',20,'FontName','Whitney-Book');
set(gca,'FontSize',20,'FontName','Whitney-Book')
legend('Modelled GPP','Measured GPP');

set(gcf,'paperPositionMode','auto') % make the print as big as the figure
print(gcf, '-dpng','-r300', ofigure1);



