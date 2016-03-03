function biochem_out = biochemical(biochem_in,Ci_input)
%
% Date: 	21 Sep 2012
% Update:   20 Feb 2013
% Update:      Aug 2013: correction of L171: Ci = Ci*1e6 ./ p .* 1E3;
%   
% Authors: 	Joe Berry and Christiaan van der Tol, contributions of others.
% Sources: 	
%           Farquhar et al. 1980, Collatz et al (1991, 1992).
%
% This function calculates:
%    - stomatal resistance of a leaf or needle (s m-1)
%    - photosynthesis of a leaf or needle (umol m-2 s-1)
%    - fluorescence of a leaf or needle (fraction of fluor. in the dark)
%
% Usage:
% biochem_out = biochemical(biochem_in)
% the function was tested for Matlab 7.2.0.232 (R2006a)
%
% Calculates net assimilation rate A, fluorescence F using biochemical model
%
% Input (units are important):
% structure 'biochem_in' with the following elements:
% Knparams   % [], [], []           parameters for empirical Kn (NPQ) model: Kn = Kno * (1+beta).*x.^alpha./(beta + x.^alpha);
%       [Kno, Kn_alpha, Kn_beta]
%  or, better, as individual fields:
%   Kno                                     Kno - the maximum Kn value ("high light")
%   Kn_alpha, Kn_beta                      alpha, beta: curvature parameters  
%
% Cs        % [ppmV or umol mol]    initial estimate of conc. of CO2 in the
%                                   ...bounary layer of the leaf
% Q         % [umol photons m-2 s-1]net radiation, PAR
% fPAR     % [0-1]                 fraction of incident light that is absorbed by the leaf (default = 1, for compatibility)
% T         % [oC or K]             leaf temperature
% eb        % [hPa = mbar]          intial estimate of the vapour pressure in leaf boundary layer
% O         % [mmol/mol]            concentration of O2 (in the boundary
%                                   ...layer, but no problem to use ambient)
% p         % [hPa]                 air pressure
% Vcmax25 (Vcmo)  % [umol/m2/s]     maximum carboxylation capacity @ 25 degC
% BallBerrySlope (m) % []           Ball-Berry coefficient 'm' for stomatal regulation
% BallBerry0 % []              (OPTIONAL) Ball-Berry intercept term 'b' (if present, an iterative solution is used)

% Type      % ['C3', 'C4']          text parameter, either 'C3' for C3 or any
%                                   ...other text for C4
% tempcor   % [0, 1]               boolean (0 or 1) whether or not
%                                   ...temperature correction to Vcmax has to be applied.
% Tparams  % [],[],[K],[K],[K]     vector of 5 temperature correction parameters, look in spreadsheet of PFTs.
%                                     Only if tempcor=1, otherwise use dummy values
% ...Or replace w/ individual values:
% slti        []              slope of cold temperature decline (C4 only)
% shti        []              slope of high temperature decline in photosynthesis
% Thl         [K]             T below which C4 photosynthesis is <= half that predicted by Q10
% Thh         [K]             T above which photosynthesis is <= half that predicted by Q10
% Trdm        [K]             T at which respiration is <= half that predicted by Q10

% effcon    [mol CO2/mol e-]  number of CO2 per electrons - typically 1/5 for C3 and 1/6 for C4

% RdPerVcmax25 (Rdparam)  % []     respiration as fraction of Vcmax25
% stressfactor [0-1]               stress factor to reduce Vcmax (for
%                                   example soil moisture, leaf age). Use 1 to "disable" (1 = no stress)
%  OPTIONAL
% Kpep25 (kp)   % [umol/m2/s]         PEPcase activity at 25 deg C (defaults to Vcmax/56
% atheta      % [0-1]                  smoothing parameter for transition between Vc and Ve (light- and carboxylation-limited photosynthesis)
% useTLforC3  % boolean              whether to enable low-temperature attenuation of Vcmax in C3 plants (its always on for C4 plants)
%
% Note: always use the prescribed units. Temperature can be either oC or K
% Note: input can be single numbers, vectors, or n-dimensional
% matrices
%
% Output:
% structure 'biochem_out' with the following elements:
% A         % [umol/m2/s]           net assimilation rate of the leaves
% Cs        % [umol/m3]             CO2 concentration in the boundary layer
% eta0      % []                    fluorescence as fraction of dark
%                                   ...adapted (fs/fo)
% rcw       % [s m-1]               stomatal resistance
% qE        % []                    non photochemical quenching
% fs        % []                    fluorescence as fraction of PAR
% Ci        % [umol/m3]             internal CO2 concentration
% Kn        % []                    rate constant for excess heat
% fo        % []                    dark adapted fluorescence (fraction of aPAR)
% fm        % []                    light saturated fluorescence (fraction of aPAR)
% qQ        % []                    photochemical quenching
% Vcmax     % [umol/m2/s]           carboxylation capacity after
%                                   ... temperature correction 

if nargin < 2
    Ci_input = [];
end
%% input
 % environmental
if isfield(biochem_in, 'Cs')
   Cs         = biochem_in.Cs;
else
    % if Cs is missing, Ci must have been supplied. Forcing Cs = Ci invalidates rcw.
   Cs         = biochem_in.Ci;
end
Q             = biochem_in.Q;
T             = biochem_in.T + 273.15*(biochem_in.T<200); % convert temperatures to K if not already
eb            = biochem_in.eb;
O             = biochem_in.O;
p             = biochem_in.p;

 % physiological
Type          = biochem_in.Type;
if isfield(biochem_in, 'Vcmax25')
    % new field names
    Vcmax25       = biochem_in.Vcmax25;
    BallBerrySlope = biochem_in.BallBerrySlope;
    RdPerVcmax25       = biochem_in.RdPerVcmax25;
else
    % old field names: Vcmo, m, Rdparam
    Vcmax25       = biochem_in.Vcmo;
    BallBerrySlope = biochem_in.m;
    RdPerVcmax25       = biochem_in.Rdparam;
end
BallBerry0 = 0.01; % default value
if isfield(biochem_in, 'BallBerry0')
    BallBerry0 = biochem_in.BallBerry0;
end
if isfield(biochem_in, 'effcon')
    effcon        = biochem_in.effcon;
elseif strcmpi('C3', Type)
    effcon =  1/5;
else
    effcon = 1/6; % C4
end

% SCOPE provides PAR as APAR, so default (for SCOPE) = 1
%    The curve-fitting GUI may not be providing APAR and should therefore explicitly set fPAR
fPAR = 1;  % fraction of incident light that is absorbed by the leaf
if isfield(biochem_in, 'fPAR')
    fPAR = biochem_in.fPAR;
end

  % physiological options
tempcor       = biochem_in.tempcor;
stressfactor  = biochem_in.stressfactor;
%model_choice  = biochem_in.Fluorescence_model;
if isfield(biochem_in, 'useTLforC3')
    useTLforC3   = biochem_in.useTLforC3;
else
    useTLforC3 = false;
end
    
  % fluoeresence
if isfield(biochem_in, 'Knparams')
    Knparams      = biochem_in.Knparams;
elseif isfield( biochem_in, 'Kn0')
    Knparams = [biochem_in.Kn0, biochem_in.Kn_alpha, biochem_in.Kn_beta];
elseif isfield(biochem_in, 'Fluorescence_model') && biochem_in.Fluorescence_model==0
    % default drought values: 
    Knparams = [5.01, 1.93, 10];
else
    % default general values (cotton dataset)
    Knparams = [2.48, 2.83, 0.114];    
end

% physiological temperature parameters: temperature sensitivities of Vcmax, etc 
if isfield(biochem_in, 'Tparams')
   Tparams = biochem_in.Tparams;
    slti        = Tparams(1);
    shti        = Tparams(2);
    Thl         = Tparams(3);
    Thh         = Tparams(4);
    Trdm        = Tparams(5);
else
    slti        = biochem_in.slti;
    shti        = biochem_in.shti;
    Thl         = biochem_in.Thl;
    Thh         = biochem_in.Thh;
    Trdm        = biochem_in.Trdm;
end

%  NOTE: kpep (kp), atheta parameters in next section

%% parameters (at optimum temperature)
Tref        = 25+273.15;        % [K]           absolute temperature at 25 oC

Kc25       = 350;              % [ubar]        kinetic coefficient (Km) for CO2 (Von Caemmerer and Furbank, 1999)
Ko25       = 450;              % [mbar]        kinetic coeeficient (Km) for  O2 (Von Caemmerer and Furbank, 1999)
spfy25     = 2600;      %  specificity (tau in Collatz e.a. 1991)
                        %     This is, in theory, Vcmax/Vomax.*Ko./Kc, but used as a separate parameter

Kpep25       = (Vcmax25/56)*1E6;      % []      (C4) PEPcase rate constant for CO2, used here: Collatz et al: Vcmax25 = 39 umol m-1 s-1; kp = 0.7 mol m-1 s-1.
if isfield(biochem_in,'Kpep25')
    Kpep25 = biochem_in.kpep;
elseif isfield(biochem_in,'kp')
    Kpep25 = biochem_in.kp;
end
if isfield(biochem_in,'atheta') 
    atheta = biochem_in.atheta;
else
    atheta      = 0.8;
end

 % electron transport and fluorescence
Kf          = 0.05;             % []            rate constant for fluorescence
%Kd          = 0.95;             % []           rate constant for thermal deactivation at Fm
Kd          = max(0.8738,  0.0301*(T-273.15)+ 0.0773);
Kp          = 4.0;              % []            rate constant for photochemisty

% note:  rhoa/Mair = L/mol (with the current units) = 24.039 L/mol
%    and  V/n = RT/P ==>  T = 292.95 K @ 1 atm (using R_hPa = 83.144621; 1 atm = 1013.25 hPa)
%   ??!!  These values are used only for rcw, however.
rhoa        = 1.2047;           % [kg m-3]       specific mass of air
Mair        = 28.96;            % [g mol-1]      molecular mass of dry air


%% convert all to bar: CO2 was supplied in ppm, O2 in permil, and pressure in mBar
ppm2bar =  1e-6 .* (p .*1E-3);
Cs          = Cs .* ppm2bar;
O           = (O * 1e-3) .* (p .*1E-3) * ~strcmp('C4',Type);    % force O to be zero for C4 vegetation (this is a trick to prevent oxygenase)
Kc25       = Kc25 * 1e-6;
Ko25       = Ko25 * 1e-3;

%% temperature corrections
qt          = 0.1 * (T-Tref) * tempcor;  % tempcorr = 0 or 1: this line dis/enables all Q10 operations
TH          = 1 + tempcor* exp(shti .* (T   -Thh));
%TH          = 1 + tempcor* exp((-220E3+703*T)./(8.314*T));
TL          = 1 + tempcor* exp(slti .* (Thl -T));

QTVc   =  2.1; % Q10 base for Vcmax and Kc
Kc          = Kc25 * 2.1.^qt;
Ko          = Ko25 * 1.2.^qt;
kpepcase          = Kpep25.* 1.8.^qt;  % "pseudo first order rate constant for PEP carboxylase WRT pi (Collatz e.a. 1992)


% jak 2014-12-04: Add TL for C3 as well, works much better with our cotton temperature dataset (A-T)
if strcmpi(Type, 'C3') && ~useTLforC3
   Vcmax = Vcmax25 .* QTVc.^qt ./TH * stressfactor;
else
   Vcmax = Vcmax25 .* QTVc.^qt ./(TL.*TH) * stressfactor;
end

% specificity (tau in Collatz e.a. 1991)
spfy        = spfy25 * 0.75 .^qt;

% "Dark" Respiration
Rd          = RdPerVcmax25 * Vcmax25 .* 1.8.^qt ./(1+exp(1.3*(T-Trdm)));

%% calculation of potential electron transport rate
po0         = Kp./(Kf+Kd+Kp);         % maximum dark photochemistry fraction, i.e. Kn = 0 (Genty et al., 1989)
Je          = 0.5*po0 .* Q .* fPAR;          % potential electron transport rate (JAK: add fPAR)

%% calculation of the intersection of enzyme and light limited curves
% this is the original Farquhar model
Gamma_star         = 0.5 .*O ./spfy; %[bar]       compensation point in absence of Rd (i.e. gamma*) [bar]

% Gamma: CO2 compensation point including Rd: solve for 0 = An = Vc - Rd, assuming Vc dominates at CO2 compensation point according to Farquar 1980. (from Leuning 1990)
if strcmp(Type, 'C3')
    Gamma = (Gamma_star .* Vcmax  +  Kc .* Rd .* (1 +O ./ Ko))  ./ (Vcmax - Rd); % C3
    MM_consts = (Kc .* (1+O./Ko));
    Vs_C3 = (Vcmax25/2) .* 1.8.^qt;
    %  minimum Ci (as fraction of Cs) for BallBerry Ci. (If Ci_input is present we need this only as a placeholder for the function call)
    minCi = 0.3;
else
    % C4
    Gamma = Vcmax - Rd;
    MM_consts = 0; % just for formality, so MM_consts is initialized
    Vs_C3 = 0;     %  the same
    minCi = 0.1;  % C4
end


%% calculation of internal CO2 concentration, photosynthesis
RH = min(1, eb./satvap(T-273.15) ); % jak: don't allow "supersaturated" air! (esp. on T curves)
warnings = [];

fcount = 0; % the number of times we called computeA()
if  ~isempty(Ci_input) 
    Ci = Ci_input; % in units of bar.
    if any(Ci_input > 1)
        % assume Ci_input is in units of ppm. Convert to bar
        Ci = Ci_input .* ppm2bar;
    end
    A =  computeA(Ci);

elseif all(BallBerry0 == 0)
    % b = 0: no need to iterate:
    Ci = BallBerry(Cs, RH, [], BallBerrySlope, BallBerry0, minCi);
    A =  computeA(Ci);
    
else
    % compute Ci using iteration (inverse-quadratic method)
    
    %  Method: "guess" Ci and then calculate the error in the guess as "Ci from B-B using A(Ci_guessed)" - Ci_guessed
    % it would be nice to use a built-in root-seeking function but fzero requires scalar inputs and outputs,
    %   and arrayfun would require converting all of the input variables to vectors (because we don't know ahead of time which inputs are scalar vs. vectors)
    % fzero uses Brent's method which requires one to (a) bracket f(x) = 0 and (b) reorder the variables on every iteration.
    %  doing it as a vector operation in MATLAB is a little tricky (see below)
    % NOTE: as best as I can determine 
    %   (1) the error function is monotonic with Ci (and often near-linear between corners caused by max() in BallBerry.m)
    %   (2) The values: x1 = Cs and x2 = the Ci calculated from A using x1, bracket error = 0;
    %      This second condition is required for Brent's method.
    tol = 1e-7;  % 0.1 ppm more-or-less

    %  Start with Ci1 = Cs and Ci2 = the first "iterated Ci" 
    Ci1 = Cs;
    [err1, Ci2] = Ci_error(Ci1, Cs, RH, minCi);  % err1 is the error for Ci1 = Cs, Ci2 ( calls computeA(Ci1)  )
    if length(Ci1) == 1 && length(err1) > 1
        % make sure Ci1 matches the size of all other Ci_n
        Ci1 = repmat(Ci1, size(err1));
    end
    err2 = Ci_error(Ci2, Cs, RH, minCi);  %  calls computeA(Ci2)
    err3 = err2;  % for the while loop; Ci3 is in case we don't enter the while loop

    % ...interlude: 
    % every test I tried in gs_tests.m confirms that the two first guesses bracket zero.
    %  this is just to further challenge the hypothesis
    %  NOTE: the algorithm does not require bracketting, though it should be much faster/well-behaved if so.
    bracketting_zero = abs(sign(err1) - sign(err2)) > 0;  % allows for err = 0
    if ~all(bracketting_zero)
        warnings = 'Not all initial Ci guesses bracket zero.';  % error?
        warning( warnings );
    end
    
    % sort the vectors so Ci1 has the lower error:
    err2_is_lower = err2 < err1;
    temp = Ci1( err2_is_lower);
    Ci1( err2_is_lower )  = Ci2( err2_is_lower);
    Ci2( err2_is_lower ) = temp; % swap in the values from Ci1
    
    temp = err1( err2_is_lower);
    err1( err2_is_lower )  = err2( err2_is_lower); % could use min here
    err2( err2_is_lower ) = temp;

    % special case: guesses that bracket gamma may encounter a sharp grief-inducing corner
    %  if we replace an endpoint with Gamma, the remaining iterations should be much simplified
    %  note: the position of the corner isn't as simple as this, but it should help some cases, at least
    bracket_gamma = ( Ci1 < Gamma & Ci2 > Gamma) | ( Ci1 > Gamma & Ci2 < Gamma);
    if any(bracket_gamma)
        % replace one endpoint with Gamma:
        Ci_Gamma = Ci2;
        if length(Gamma) == 1
            Gamma = repmat(Gamma, size(Ci1));
        end
        Ci_Gamma(bracket_gamma) = Gamma(bracket_gamma);
        errGamma = Ci_error(Ci_Gamma, Cs, RH, minCi);  %calls computeA(Ci_Gamma)
       
        % sort the result into Ci1, Ci2
        errGamma_neg = bracket_gamma & errGamma < 0;
        err1(errGamma_neg) = errGamma(errGamma_neg);
        Ci1(errGamma_neg) = Gamma(errGamma_neg);
        
        errGamma_pos = bracket_gamma & errGamma >= 0;
        err2(errGamma_pos) = errGamma(errGamma_pos);
        Ci2(errGamma_pos) = Gamma(errGamma_pos);
        
        % note: I didn't update err3, but that should be fine (there's probably huge error on everything)
    end
    
    % at this point err1, Ci1 correspond to the position with error < 0
    %      and      err2, Ci2 correspond with error >= 0
    
    % initialize Ci3_sorted so we don't reallocate on every iteration. (not sure it does anything)
    Ci3_sorted = zeros(size(Ci2));
    Ci3_sorted(1) = -1;
    err3_sorted = zeros(size(Ci2));
    Ci3 = [];  % signal first iteration
    while any( abs(err3) > tol )
        
        if ~isempty(Ci3)  %not first iteration, swap values in/out of Ci3_sorted
            % sort the result into the existing vectors (assume we're bracketting zero so Ci1 < 0 and Ci2 > 0)
            % note that values in Ci3 are, by design, closer to zero than the values in Ci1 or Ci2 with the same sign.
            % and, more generally, Ci1 < Ci3 < Ci2
            % Ci3_sorted, however, will be either < Ci1 or > Ci2 when we're done
            Ci3_is_neg = sign(err3) < 0;
            % swap negative values with Ci1
            Ci3_sorted( Ci3_is_neg ) = Ci1( Ci3_is_neg );
            err3_sorted( Ci3_is_neg ) = err1( Ci3_is_neg );
            Ci1( Ci3_is_neg ) = Ci3( Ci3_is_neg );
            err1( Ci3_is_neg ) = err3( Ci3_is_neg );

            % swap positive values with Ci2
            Ci3_isnt_neg = ~Ci3_is_neg;  % this is repeated often enough to warrant its own variable (saves ~4% overall )
            Ci3_sorted( Ci3_isnt_neg ) = Ci2( Ci3_isnt_neg );
            err3_sorted( Ci3_isnt_neg ) = err2( Ci3_isnt_neg );
            Ci2( Ci3_isnt_neg ) = Ci3( Ci3_isnt_neg );
            err2( Ci3_isnt_neg ) = err3( Ci3_isnt_neg );
        end
        
        if Ci3_sorted(1) == -1  %first iteration   (was: isempty(Ci3_sorted)  )
            % compute the next point based on the last two using the secant method
                %      linear solution for y = 0 (x-intercept, given two points):
                %      (y - y1) = m(x - x1) + b:  m = (y2-y1)/(x2-x1);  b = y1 - m x1 = y2 - m x2;  
                %      0 = mx + b => x = -b/m  = -(y1 - m x1)/m =  -y1/m + x1 = x1 - y1 (x2 - x1)/(y2 - y1)
                % note: it's possible to get negative Ci values from this solution (though not if we're bracketting zero!)
            Ci3 = Ci1 - err1 .* (Ci2 - Ci1) ./ (err2 - err1);  % or Ci2 - err2*(Ci2 - Ci1)/(err2 - err1);
            y_underflow = err2 == err1;
            if any( y_underflow )
                assert(all(Ci2( y_underflow ) == Ci1( y_underflow )), 'Two distinct values of Ci have the same error value while iterating for g_s.');
                % if err2 == err1, Ci1 should equal Ci2, so we're done (with that Ci)
                Ci3( y_underflow ) = Ci2( y_underflow ) ;
            end
%          use_mean = err2 - err1 > 1e-4 & (-err1 > 2*err2 | err2 > -err1*2);  % recall we've arranged: err1 < 0, err2 > 0;
%          Ci3_mean = (Ci2(use_mean) + Ci1(use_mean))./2;
%          Ci3(use_mean) = Ci3_mean;
        else
            % three values are available, use inverse quadratic method
            
            % compute the inverse quadratic estimate (note: Ci3_sorted is outside of (Ci1, Ci2) )
            Ci3 = Ci1.*err2.*err3_sorted ./ ((err1 - err2).*(err1 - err3_sorted)) ...
                  + Ci2.*err1.*err3_sorted ./ ((err2 - err1).*(err2 - err3_sorted))  ...
                  + Ci3_sorted.*err1.*err2 ./ ((err3_sorted - err1).*(err3_sorted - err2));
            %y_underflow = (err2 == err1); % this should always be true if the following are true: | (err3_sorted == err1) | (err3_sorted == err2);
            y_underflow = isnan(Ci3) | ~isfinite(Ci3);
            if any( y_underflow )
                %assert(all(Ci2( y_underflow ) == Ci1( y_underflow )), 'Two distinct values of Ci have the same error value while iterating for g_s.');
                % if err2 == err1, Ci1 should equal Ci2, so we're done (with that Ci)
                err1_is_err2 = err2 == err1;
                use_Ci2 = err1_is_err2 | err3_sorted == err2;
                Ci3( use_Ci2 ) = Ci2( use_Ci2 ) ;
                % now only use Ci1 if Ci1 ~= Ci2:
                use_Ci1 = ~err1_is_err2 & err3_sorted == err1;
                Ci3( use_Ci1 ) = Ci1( use_Ci1 ) ;
                assert( all(  ~isnan(Ci3) & isfinite(Ci3) ), 'NaN or Inf in Ci3 were not fully removed.');
            end
        end
        Ci3 = max(0, Ci3);
        % A problem with the secant method is that it sticks to one side of the zero line, depending on the concavity of the function
        % So one side creeps towards zero while the other is stationary.  We can converge more quickly by trying to balance
        % the adjustments on both sides (thus reducing distance from zero, and therefore curvature, for both points).
        % Once we're pretty close, however, the curvature is likely small, and it will be better to use the second.
        %   I've provisionally defined 'small' as < 100 ppm apart (alt, define small relative to zero?)
%         use_mean = err2 - err1 > 1e-4 & (-err1 > 2*err2 | err2 > -err1*2);  % recall we've arranged: err1 < 0, err2 > 0;
%         Ci3_mean = (Ci2(use_mean) + Ci1(use_mean))./2;
%         Ci3(use_mean) = Ci3_mean;
        
        err3 = Ci_error(Ci3, Cs, RH, minCi);  % calls  computeA(Ci3)
    end
    if isempty(Ci3)
        % we never entered the while loop, so Ci2 is the correct result (and is consistent with the last call to computeA() ).
        Ci = Ci2;
    else
        % use Ci3 from the while loop:
        Ci = Ci3;
    end
    % A has been computed in the last call to Ci_error
end

%%
    function [err, Ci_out] = Ci_error(Ci_in, Cs, RH, minCi)
        % compute the difference between "guessed" Ci (Ci_in) and Ci computed using BB after computing A
        A = computeA(Ci_in);
        A_bar = A .* ppm2bar;
        [Ci_out, ~] = BallBerry(Cs, RH, A_bar, BallBerrySlope, BallBerry0, minCi, Ci_input); %[Ci_out, gs]
        err = Ci_out - Ci_in;
    end

%%
    function [A, biochem_out] = computeA(Ci)
        % global: Type, Vcmax, Gamma_star, MM_consts, Vs_C3, effcon, Je, atheta, Rd    %Kc, O, Ko, Vcmax25, qt
        if strcmpi('C3', Type)
            %[Ci, gs] = BallBerry(Cs, RH, A_bar, BallBerrySlope, BallBerry0, 0.3, Ci_input);
            Vc          = Vcmax.*(Ci-Gamma_star)./(MM_consts + Ci);  % MM_consts = (Kc .* (1+O./Ko)) % doesn't change on iteration.
            Vs          = Vs_C3; % = (Vcmax25/2) .* 1.8.^qt;    % doesn't change on iteration.
           %effcon      = 0.2;
            CO2_per_electron = (Ci-Gamma_star)./(Ci+2*Gamma_star) .* effcon;

        else  %C4
            %[Ci, gs] = BallBerry(Cs, RH, A_bar, BallBerrySlope, BallBerry0, 0.1, Ci_input);
            Vc          = Vcmax;
            Vs          = kpepcase.*Ci;
            %effcon      = 0.17;                    % Berry and Farquhar (1978): 1/0.167 = 6
            CO2_per_electron = effcon; % note: (Ci-Gamma_star)./(Ci+2*Gamma_star) = 1 for C4 (since O = 0); this line avoids 0/0 when Ci = 0
        end
        Ve          = Je .* CO2_per_electron;

        % find the smoothed minimum of Ve, Vc = V, then V, Vs
%         [a1,a2]     = abc(atheta,-(Vc+Ve),Vc.*Ve);
%         % select the min or max  depending on the side of the CO2 compensation point
%         %  note that Vc, Ve < 0 when Ci < gam, so max will select the root that is closer to 0
%         V           = min(a1,a2).*(Ci>Gamma_star) + max(a1,a2).*(Ci<=Gamma_star);
%         [a1,a2]     = abc(0.98,-(V+Vs),V.*Vs);
%         Ag          = min(a1,a2);
        V           = min_root(atheta,-(Vc+Ve),Vc.*Ve);
        Ag          = min_root(0.98,-(V+Vs),V.*Vs);
        A           = Ag - Rd;
        %A           = Ag; %XY test for gross A
        if nargout > 1
            biochem_out.A = A;
            biochem_out.Ag = Ag;
            biochem_out.Vc = Vc;
            biochem_out.Vs = Vs;
            biochem_out.Ve = Ve;
            biochem_out.CO2_per_electron = CO2_per_electron;
        end
        fcount = fcount + 1; % # of times we called computeA

    end

% (ppm2bar), A_bar
%tic;
%toc
%fprintf('Ball-Berry converged in %d steps (largest_diff = %.4g)\n', counter, largest_diff/ mean(ppm2bar));
%% Compute A, etc.

% note: the following sets a bunch of "global" values in the nested function. Prob better to use [A biochem_out] = ....
%A =  computeA(Ci);  % done above
[~, gs] = BallBerry(Cs, RH, A .* ppm2bar, BallBerrySlope, BallBerry0, minCi, Ci);

Ja          = Ag ./ CO2_per_electron;        % actual electron transport rate

 % stomatal resistance 
%old: rcw         = 0.625*(Cs-Ci)./A *rhoa/Mair*1E3  ./ ppm2bar; %  * 1e6 ./ p .* 1E3;
% if BallBerry0 == 0  %if B-B intercept was specified, then we computed gs "correctly" above and don't need this.
%     rcw(A<=0 & rcw~=0)   = 0.625*1E6;
% end
%rcw         = (1./gs) *rhoa/Mair*1E3  ./ ppm2bar; %  * 1e6 ./ p .* 1E3;
rcw      =  (rhoa./(Mair*1E-3))./gs;
%% fluorescence (Replace this part by Magnani or other model if needed)
ps          = po0.*Ja./Je;               % this is the photochemical yield
nanPs = isnan(ps);
if any(nanPs)
    if numel(po0) == 1
        ps(nanPs) = po0;
    else
        ps(nanPs) = po0(nanPs);  % happens when Q = 0, so ps = po0 (other cases of NaN have been resolved) 
    end
end
ps_rel   = max(0,  1-ps./po0);       % degree of light saturation: 'x' (van der Tol e.a. 2014)

[eta,qE,qQ,fs,fo,fm,fo0,fm0,Kn]    = Fluorescencemodel(ps, ps_rel, Kp,Kf,Kd,Knparams);
Kpa         = ps./fs*Kf;

%% convert back to ppm
Ci          = Ci*1e6 ./ p .* 1E3;
%Cs          = Cs*1e6 ./ p .* 1E3;

%% Collect outputs

biochem_out.A       = A;
biochem_out.Ag      = Ag; %XY added Gross Photosynthesis as the output
biochem_out.Ci      = Ci;
biochem_out.rcw     = rcw;
biochem_out.gs      =  gs;
%  this would be the same if we apply the rcw(A<=0) cutoff:
%biochem_out.gs      = BallBerrySlope.*A.*RH./Cs; % mol/m2/s no intercept term.
biochem_out.RH      =  RH;
biochem_out.warnings = warnings;
biochem_out.fcount = fcount;  % the number of times we called computeA()
%fprintf('fcount = %d\n', fcount);

biochem_out.Vcmax   = Vcmax;
biochem_out.Vc = Vc;  % export the components of A for diagnostic charts
biochem_out.Ve = Ve;
biochem_out.Vs = Vs;
biochem_out.Rd = Rd;

biochem_out.Ja      = Ja;
biochem_out.ps      = ps; % photochemical yield
biochem_out.ps_rel  = ps_rel;   % degree of light saturation 'x' (van der Tol e.a. 2014)

 % fluoresence outputs:
biochem_out.Kd      = Kd;  % K_dark(T)
biochem_out.Kn      = Kn;  % K_n(T)
                           % Kf = 0.05 (const)
                           % Kp = 4.0 (const) ??
biochem_out.Kp      = Kpa;
biochem_out.eta     = eta;
biochem_out.qE      = qE;
biochem_out.fs      = fs;  % keep this for compatibility with SCOPE
biochem_out.ft      = fs;  % keep this for the GUI ft is a synonym for what we're calling fs
biochem_out.SIF     = fs .* Q;
biochem_out.fo0     = fo0;
biochem_out.fm0     = fm0;
biochem_out.fo      = fo;
biochem_out.fm      = fm;
biochem_out.qQ      = qQ;
return;

end  % end of function biochemical
 
%% quadratic formula
function [x2,x1] = abc(a,b,c)
    if a == 0
        x1      = -c./b;
        x2      = x1;
    else
        det_root = sqrt(b.^2 - 4.*a.*c);
        x1      = (-b+det_root)./(2.*a);
        x2      = (-b-det_root)./(2.*a);
    end
end %of abc formula

%% quadratic formula, root of least magnitude
function x = min_root(a,b,c)
    if a == 0  % note: this works because 'a' is a scalar parameter!
        x      = -c./b;
    else
        %disc_root = sqrt(b.^2 - 4.*a.*c); % square root of the discriminant (doesn't need a separate line anymore)
        %  in MATLAB (2013b) assigning the intermediate variable actually slows down the code! (~25%)
        x = (-b + sign(b).* sqrt(b.^2 - 4.*a.*c))./(2.*a);
    end
end %of min_root of quadratic formula


%% Ball Berry Model
function [Ci, gs] = BallBerry(Cs, RH, A, BallBerrySlope, BallBerry0, minCi, Ci_input)
%  Cs  : CO2 at leaf surface
%  RH  : relative humidity
%  A   : Net assimilation in 'same units of CO2 as Cs'/m2/s
% BallBerrySlope, BallBerry0, 
% minCi : minimum Ci as a fraction of Cs (in case RH is very low?)
% Ci_input : will calculate gs if A is specified.
if nargin > 6 && ~isempty(Ci_input)
    % Ci is given: try and compute gs
    Ci = Ci_input;
    gs = [];
    if ~isempty(A)
        gs = gsFun(Cs, RH, A, BallBerrySlope, BallBerry0);
    end
elseif all(BallBerry0 == 0) || isempty(A)
    % EXPLANATION:   *at equilibrium* CO2_in = CO2_out => A = gs(Cs - Ci) [1]
    %  so Ci = Cs - A/gs (at equilibrium)                                 [2]
    %  Ball-Berry suggest: gs = m (A RH)/Cs + b   (also at equilib., see Leuning 1990)
    %  if b = 0 we can rearrange B-B for the second term in [2]:  A/gs = Cs/(m RH)
    %  Substituting into [2]
    %  Ci = Cs - Cs/(m RH) = Cs ( 1- 1/(m RH)  [ the 1.6 converts from CO2- to H2O-diffusion ]
    Ci      = max(minCi*Cs,  Cs.*(1-1.6./(BallBerrySlope .* RH)));
    gs = [];
else
    %  if b > 0  Ci = Cs( 1 - 1/(m RH + b Cs/A) )
    % if we use Leuning 1990, Ci = Cs - (Cs - Gamma)/(m RH + b(Cs - Gamma)/A)  [see def of Gamma, above]
    % note: the original B-B units are A: umol/m2/s, ci ppm (umol/mol), RH (unitless)
    %   Cs input was ppm but was multiplied by ppm2bar above, so multiply A by ppm2bar to put them on the same scale.
    %  don't let gs go below its minimum value (i.e. when A goes negative)
    gs = gsFun(Cs, RH, A, BallBerrySlope, BallBerry0);
    Ci = max(minCi*Cs,  Cs - 1.6 * A./gs) ; % mean of the new and the old
    %Ci = max(.3*Cs,  Cs - 1.6./(BallBerrySlope.*RH ./Cs  + BallBerry0 ./ ( A*1e-6  .*  p*1E-3 )) );
end

end % function

function gs = gsFun(Cs, RH, A, BallBerrySlope, BallBerry0)
% add in a bit just to avoid div zero. 1 ppm = 1e-6
gs = max(BallBerry0,  BallBerrySlope.* A .* RH ./(Cs+10e-9)  + BallBerry0);
% alt:
%     % eliminate infinities
%     zeroCs = Cs == 0;
%     if any(zeroCs) 
%         gs( zeroCs ) = repmat(BallBerry0, sum(zeroCs), 1);
%     end

end


%% Fluorescence model
function [eta,qE,qQ,fs,fo,fm,fo0,fm0,Kn] = Fluorescencemodel(ps,x, Kp,Kf,Kd,Knparams)
  % note: x isn't strictly needed as an input parameter but it avoids code-duplication (of po0) and it's inherent risks.
  
    Kno = Knparams(1);
    alpha = Knparams(2);
    beta = Knparams(3);

    % switch model_choice
    %     case 0, % drought 
    %         Kno = 5.01;
    %         alpha = 1.93;
    %         beta = 10;
    %         %Kn          = (6.2473 * x - 0.5944).*x; % empirical fit to Flexas' data
    %         %Kn          = (3.9867 * x - 1.0589).*x;  % empirical fit to Flexas, Daumard, Rascher, Berry data
    %     case 1, healthy (cotton)
    %         Kno = 2.48;
    %         alpha = 2.83;
    %         beta = 0.114;
    %         %p = [4.5531;8.5595;1.8510];
    %         %Kn   = p(1)./(p(3)+exp(-p(2)*(x-.5)));
    % end

    % using exp(-beta) expands the interesting region between 0-1
    %beta = exp(-beta);
    x_alpha = x.^alpha; % this is the most expensive operation in this fn; doing it twice almost doubles the time spent here (MATLAB 2013b doesn't optimize the duplicate code)
    Kn = Kno * (1+beta).* x_alpha./(beta + x_alpha);

    %Kn          = Kn .* Kd/0.8738;          % temperature correction of Kn similar to that of Kd

    fo0         = Kf./(Kf+Kp+Kd);           % dark adapted fluorescence yield Fo -- this should be the quantum yield of minimal Chl a fluorescence
    fo          = Kf./(Kf+Kp+Kd+Kn);           % XY: light adapted fluorescence yield Fo'
    fm          = Kf./(Kf   +Kd+Kn);        % light adapted fluorescence yield Fm' (XY: Eq.13 in vdT 2014)
    fm0         = Kf./(Kf   +Kd);        %  XY: dark-adapted fluorescence yield Fm
    fs          = fm.*(1-ps);            %(Eq.12 in vdT 2014)
    eta         = fs./fo0;
    qQ          = 1-(fs-fo)./(fm-fo);       % photochemical quenching
    qE          = 1-(fm-fo)./(fm0-fo0);     % non-photochemical quenching

    eta         = eta*(1+5)/5 - 1/5;

end


