% This script generates parts of the results from:
% "The uncertainty of the experimentally-measured momentum coefficient:
% Guidlines on how to accurately estimate Cmu," by Richard Semaan
%
%
% Copyright (c) 2020, Richard Semaan
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution.
% 
% * Neither the name of SCOUT nor the names of its
%   contributors may be used to endorse or promote products derived from
%   this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% Author: Richard Semaan
% Affiliation: Technische Universität Braunschweig
% Contact: r.semaan@tu-bs.de
% Release date: Apr. 18, 2020
% Version: v 1.0
% 
% Installation
% -------------------------
% - No installation needed
% 
% Input
% -------------------------
% - Selection among the various approaches and 3 scenarios
% - Individual random and systematic uncertainties of various variables
% - Gas constant, pressure ratio, jet exit pressure, blowing intensity*, nominal height* (* denotes variables that are only needed for the UPCs)
% 
% Output
% -------------------------
% - Uncertainty of Cmu
% - UPC of all deopendent variables
% 
% ===============================================================================================
clear; close all; clc

%% Input
compressible = 'no';   % options: 'yes' or 'no'
method = 'pressure';    % options: 'combined', 'pressure', or 'massflow'
scenario = 3;           % options: 1, 2, or 3

%% Constants
g = 1.4;        % Gas constant
PR = 1.5;       % pressure ratio
Pj = 100000;    % Exact value is not critica for the current analysis
Pp = PR.*Pj;
cmu = 0.035;    % mid range blowing intensity
h = 0.04;       % nominal height
%% UMFs
switch compressible
    case 'yes'
        switch method
            case 'combined'
                UMFm = 1;
                UMFTp = 0.5;
                UMFTj = 0;
                UMFPj = (g-1).*Pj./(2.*g.*(Pj-Pp.*(Pj./Pp).^(1./g))) ;
                UMFPp = -(g-1).*Pj./(2.*g.*(Pj-Pp.*(Pj./Pp).^(1./g))) ;
                UMFPi = -1;
                UMFh = 0;
                dcmudh = 0;
                dcmudPj = 1; % Not 1, but does not matter since dcmudh=0
            case 'pressure'
                UMFm = 0;
                UMFTp = 1;
                UMFTj = -1;
                UMFPj = ( (2.*g-1).*Pj - Pp.*g.*(Pj./Pp).^(1./g) )./(g.*(Pj-Pp.*(Pj./Pp).^(1./g))) ;
                UMFPp = - (g-1).*Pj./(g.*(Pj-Pp.*(Pj./Pp).^(1./g))) ;
                UMFPi = -1;
                UMFh = 1;
                dcmudh = cmu/h;
                dcmudPj = UMFPj*cmu/Pj; %-cmu*( (gPj^2-g*Pj*Pp*(Pj/Pp)^(1/g))/(Pj-2*g*Pj*Pp*(Pj/Ppl)^(1/g)) ); 
            case 'massflow'
                UMFm = 2;
                UMFTp = 0;
                UMFTj = 1;
                UMFPj = -1 ;
                UMFPp = 0 ;
                UMFPi = -1;
                UMFh = -1; 
                dcmudh = -cmu/h;
                dcmudPj = -cmu/Pj;
            otherwise
                error('Invalid method input.')
        end
    case 'no'
        switch method
            case 'combined'
                UMFm = 1;
                UMFTp = 0;
                UMFTj = 0;
                UMFPj = - Pj./(2.*(Pp-Pj)) ;
                UMFPp =  Pp./(2.*(Pp-Pj)) ;
                UMFPi = -1;
                UMFh = 0;
                dcmudh = 0;
                dcmudPj = 1; % Not 1, but does not matter since dcmudh=0
            case 'pressure'
                UMFm = 0;
                UMFTp = 0;
                UMFTj = 0;
                UMFPj = - Pj./((Pp-Pj)) ;
                UMFPp =  Pp./((Pp-Pj)) ;
                UMFPi = -1;
                UMFh = 1; 
                dcmudh = cmu/h;
                dcmudPj = UMFPj*cmu/Pj; 
            case 'massflow'
                UMFm = 2;
                UMFTp = 0;
                UMFTj = 0;
                UMFPj =  0;
                UMFPp =  0;
                UMFPi = -1;
                UMFh = -1; 
                dcmudh = 1;     % Not 1, but does not matter since dcmudPj=0
                dcmudPj = 0;
            otherwise
                error('Invalid method input.')
        end  
    otherwise
            error('Invalid compressibility input.')
end

%% Uncertainties
switch scenario
    case 1
        % m
        sm = 0.1/100;           % Random uncertainty of m
        bDm = 1/100;            % Device systematic uncertainty of mass flowmeter
        bSm = 0;                % Setup systematic uncertainty of the mass flowmeter
        % T
        sT = 0.02/100;          % Random uncertainty of the temperature sensors
        bDT = 0.75/100;         % Device systematic uncertainty of the thermocouples
        bST = 0;                % Setup systematic uncertainty of the thermocouples
        % Pp
        sPp = 0.09/100;         % Random uncertainty of the plenum pressure sensor 
        bDPp = 0.05/100;        % Device systematic uncertainty of Ppl
        bSPp = 0;               
        % Pi
        sPi = 0.01/100;         % Random uncertainty of the freestream pressure
        bDPi = 0.02/100;        % Device systematic uncertainty of freestream pressure
        bSPi = 0;               % Setup systematic uncertainty of the freestream pressure
        % h
        sh = 0;                 % Random uncertainty of the lip height distribution
        bDh = 0;                % Device systematic uncertainty of the lip height
        bSh1 = 10/100;          % First setup systematic uncertainty of h
        bSh2 = 5/100;           % Second setup systematic uncertainty of h
        bSh = sqrt(bSh1^2 + bSh2^2);
        % Pj
        sPj = 0.02/100;         % Random uncertainty of the jet exit pressure
        bDPj = 0.05/100;        % Device systematic uncertainty of Pj
        bSPj1 = 10/100;         % First setup systematic uncertainty Pj
        bSPj2 = bSh;            % Second setup systematic uncertainty of Pj = lip heigh uncertainty
        bSPj = sqrt(bSPj1^2 + bSPj2^2);
        % h-Pj
        bhPj = bSPj2*bSh1 + bSPj2*bSh2; % Correlated uncertainty
    case 2
        % m
        sm = 0.1/100;
        bDm = 1/100;
        bSm = 0;
        % T
        sT = 0.02/100;
        bDT = 0.75/100;
        bST = 0;   
        % Pp
        sPp = 0.09/100;
        bDPp = 0.05/100;
        bSPp = 0;  
        % Pi
        sPi = 0.01/100;
        bDPi = 0.02/100;
        bSPi = 0;  
        % h
        sh = 0;
        bDh = 0;
        bSh1 = 10/100;
        bSh2 = 5/100;
        bSh = sqrt(bSh1^2 + bSh2^2);
        % Pj
        sPj = 0.02/100;
        bDPj = 0.05/100;
        bSPj1 = 0.05/100; 
        bSPj2 = bSh;
        bSPj = sqrt(bSPj1^2 + bSPj2^2);
        % h-Pj
        bhPj = bSPj2*bSh1 + bSPj2*bSh2;        
    case 3
        % m
        sm = 0.1/100;
        bDm = 1/100;
        bSm = 0;
        % T
        sT = 0.02/100;
        bDT = 0.75/100;
        bST = 0;    
        % Pp
        sPp = 0.09/100;
        bDPp = 0.05/100;
        bSPp = 0;  
        % Pi
        sPi = 0.01/100;
        bDPi = 0.02/100;
        bSPi = 0;  
        % h
        sh = 0;
        bDh = 0;
        bSh1 = 5/100;
        bSh2 = 2.5/100;
        bSh = sqrt(bSh1^2 + bSh2^2);
        % Pj
        sPj = 0.02/100;
        bDPj = 0.05/100;
        bSPj1 = 0.05/100; 
        bSPj2 = bSh;
        bSPj = sqrt(bSPj1^2 + bSPj2^2);
        % h-Pj
        bhPj = bSPj2*bSh1 + bSPj2*bSh2;        
    otherwise
        error('Invalid scenario input.')
end

%% Combined systematic uncertainties
bm = sqrt(bDm^2 + bSm^2);
bT = sqrt(bDT^2 + bST^2);
bPj = sqrt(bDPj^2 + bSPj^2);
bPp = sqrt(bDPp^2 + bSPp^2);
bPi = sqrt(bDPi^2 + bSPi^2);
bh = sqrt(bDh^2 + bSh^2);
%% Combined uncertainty
um = sqrt(sm^2 + bm^2);
uT = sqrt(sT^2 + bT^2);
uPj = sqrt(sPj^2 + bPj^2);
uPp = sqrt(sPp^2 + bPp^2);
uPi = sqrt(sPi^2 + bPi^2);
uh = sqrt(sh^2 + bh^2);
%% Uncertainty propagation
% random uncertainty
s = sqrt( UMFm^2*sm^2   + UMFTj^2*sT^2   + UMFTp^2*sT^2 + UMFPj^2*sPj^2 ...
  + UMFPp^2*sPp^2 + UMFPi^2*sPi^2  + UMFh^2*sh^2 );
% systematic uncertainty
b = sqrt( UMFm^2*bm^2   + UMFTj^2*bT^2   + UMFTp^2*bT^2 + UMFPj^2*bPj^2 ...
  + UMFPp^2*bPp^2 + UMFPi^2*bPi^2  + UMFh^2*bh^2 ...
  + 2*dcmudh*dcmudPj*bhPj );
% Combined uncertainty
u = sqrt(s^2 + b^2);
% Combined expanded uncertainty
U = 2*u*100;
sprintf('The total relative expanded uncertainty on Cmu is: %f',round(U,2))

%% UPCs
UPCm  = 100*UMFm^2*um^2/u^2;
UPCTj = 100*UMFTj^2*uT^2/u^2;
UPCTp = 100*UMFTp^2*uT^2/u^2;
UPCPj = 100*UMFPj^2*uPj^2/u^2;
UPCPp = 100*UMFPp^2*uPp^2/u^2;
UPCPi = 100*UMFPi^2*uPi^2/u^2;
UPCh  = 100*UMFh^2*uh^2/u^2;
UPChPj = 100 - (UPCm + UPCTj + UPCTp + UPCPj + UPCPp + UPCPi + UPCh);

sprintf('The UPCs are:')
sprintf('UPC for m: %f',round(UPCm,2))
sprintf('UPC for Tj: %f',round(UPCTj,2))
sprintf('UPC for Tp: %f',round(UPCTp,2))
sprintf('UPC for Pj: %f',round(UPCPj,2))
sprintf('UPC for Pp: %f',round(UPCPp,2))
sprintf('UPC for Pi: %f',round(UPCPi,2))
sprintf('UPC for h: %f',round(UPCh,2))
sprintf('UPC for correlated uncertainty of h and Pj: %f',round(UPChPj,2))
