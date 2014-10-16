% %% Sensitivity analysis of SCOPE model
% %  SCOPE version: v1.52
% 
% clc
% clear all
% 
% sens_analysis = 1;
% 
% %  1.generate random parameters as inputs to SCOPE using Latin Hypercube
% 
% p   = 300;   % Number of runs
% N   = 3;   % Number of dimensions
% lb  = [1 1 1]; % lower bounds 100 30
% ub  = [60 100 5]; % upper bounds 800 70
% X   = lhsdesign(p,N,'criterion','correlation');  %Latin Hypercube Sampling
% D   = bsxfun(@plus,lb,bsxfun(@times,X,(ub-lb))); %Generate parameters
% 
% %   2.read input text file and add the generated parameters into the
% %   correct locations
% X                           = textread('../inputdata.txt','%s');
% N                           = str2double(X);
% V1                          = assignvarnames;
% [V1.Val]                    = deal(-999);
% 
% 
% parameter_strings = {'Cab';'Vcmo';'LAI'};%;'Rin';'tts'};
% kk = zeros(numel(parameter_strings),1);
% 
% for ii = 1:numel(parameter_strings)
%     k           = find(strcmp({V1.Name},parameter_strings(ii)));
%     V1(k).Val   = D(:,ii);
%     kk(ii)      = k;
% end
% 
% 
% 
% %   3. Run SCOPE_mac_linux
% %   I made some modifications so that it can take values from this program
% 
% SCOPE_mac_linux


%   4. Extract the output values, e.g., photosynthesis, water flux etc. for
%   Random Forest analysis

directories         =   dir('../output/*-*');                               % [XY] I made changes to this for mac system, as there are ., .., and .DS_store
[time_value_s,I]    =   sort([directories(:).datenum]);                     % [XY] Changes are made for mac for lines 20-22
Directory           =   directories(I(end)).name;                           % 

wl                  = dlmread(['../output/' Directory '/wl.dat'],'',2,0);
fEnergy             = dlmread(['../output/' Directory '/fluorescence.dat'],'',2,0);
D                   = dlmread(['../output/' Directory '/pars_and_input_short.dat'],'',1,0);


wl                  = wl(1,wl(1,:) >=640 & wl(1,:)<=850);
f760                = fEnergy(:,wl(1,:) == 760);

clearvars -except 'D' 'f760'

% Normalize D

% normD  = max(D,[],1)-min(D,[],1);
% normD1 = repmat(normD, [length(D) 1]);
% minD   = repmat(min(D,[],1), [length(D) 1]);
% normD2 = (D-minD)./normD1;



% sub_tr = 1:floor(1/2*numel(f760));
% sub_te = floor(1/2*numel(f760)):numel(f760);
extra_options.importance = 1;
model1 = regRF_train(D,f760,500,5,extra_options)
%model  = regRF_train(D,f760,extra_options)

%[test_set_predictions,chain_stats,MODEL,mxsx,mnmx,SIG2] = bayes_lm([f760(sub_tr),D(sub_tr,:)],[f760(sub_te),D(sub_te,:)]);
% extra_options.importance = 1;
% model = classRF_train(D,f760,extra_options)2