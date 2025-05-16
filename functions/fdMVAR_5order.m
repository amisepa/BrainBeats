%% FREQUENCY DOMAIN MVAR ANALYSIS
% 
% 
% INPUTS:
%   Am = [A(1)...A(p)]: M*pM matrix of the MVAR model coefficients (strictly causal model)
%   Su = M*M covariance matrix of the input noises
%   N  = number of points for calculation of the spectral functions (nfft)
%   Fs = sampling frequency
% 
% OUTPUTS:
%   DC  = Directed Coherence (Eq. 11). Directional influence from one signal
%       to another, normalized by total input to the target.
%   DTF = Directed Transfer Function (Eq. 11 but with sigma_i=sigma_j for  each i,j):
%       Directional flow from source to target, based on the transfer function, 
%       assuming equal noise across signals.
%   PDC = Partial Directed Coherence (Eq. 15 but with sigma_i=sigma_j for each i,j)
%       Direct influence from one signal to another, normalized by total 
%       output from the source, with equal noise.
%   GPDC = Generalized Partial Directed Coherence (Eq. 15): Like PDC but 
%       accounts for actual noise levels in each signal.
%   COH = Coherence (Eq. 3): Strength of connection between two signals at 
%       each frequency.
%   PCOH = Partial Coherence (Eq. 3): Connection between two signals after
%       after removing the influence of all other signals, isolating unique shared variance.
%   H   = Tranfer Function Matrix (Eq. 6): System’s frequency response from
%       inputs to outputs.
%   S   = Spectral Matrix (Eq. 7): Power and shared power between signals
%       across frequencies.
%   P   = Inverse Spectral Matrix (Eq. 7): Inverse of the spectral matrix, 
%       shows conditional dependencies.
%   f   = frequency vector
% 
% Copyright (C), Cedric Cannard, BrainBeats 2024
% 
% PLEASE CITE THE FOLLOWING REFERENCE WHEN USING THIS CODE:
%   Faes & Nollo (2011). Multivariate Frequency Domain Analysis of Causal Interactions in Physiological Time Series, Biomedical Engineering, Trends in Electronics, Communications and Software

function [DC,DTF,PDC,GPDC,IPDC,COH,PCOH,PCOH2,H,S,P,f] = fdMVAR_5order(Am,Su,N,Fs)

M = size(Am,1); % Am has dim M*pM
p = 5; % p is the order of the MVAR model

if nargin<2, Su = eye(M,M); end % if not specified, we assume uncorrelated noises with unit variance as inputs
if nargin<3, N = 512; end
if nargin<4, Fs= 1; end
if all(size(N)==1)	 %if N is scalar
    f = (0:N-1)*(Fs/(2*N)); % frequency axis
else            % if N is a vector, we assume that it is the vector of the frequencies
    f = N; N = length(N);
end

% s = exp(1i*2*pi*f/Fs); % vector of complex exponentials
z = 1i*2*pi/Fs;


%% Initializations: spectral matrices have M rows, M columns and are calculated at each of the N frequencies
H=zeros(M,M,N); % Transfer Matrix
S=zeros(M,M,N); % Spectral Matrix
P=zeros(M,M,N); % Inverse Spectral Matrix
COH=zeros(M,M,N); %Coherence
PCOH=zeros(M,M,N); %Partial Coherence - defined as Dahlhaus 2000
PCOH2=zeros(M,M,N); %Partial Coherence - defined as Yacoub 1970
DC=zeros(M,M,N); % directed coherence - defined as Baccala 1998
DTF=zeros(M,M,N); % directed transfer function - Defined as Kaminski 1991
PDC=zeros(M,M,N); % PDC Baccala Sameshima 2001
GPDC=zeros(M,M,N); %generalized PDC Baccalà Cardiff2007
IPDC=zeros(M,M,N);%isolated coherence
tmp1=zeros(M,1); %denominator for DC (column!)
tmp2=tmp1; %denominator for DC (column!)
tmp3=tmp1'; tmp4=tmp1'; %denominators for PDC (row!)

A = [eye(M) -Am]; % matrix from which M*M blocks are selected to calculate spectral functions
% invSu = inv(Su);

% Define the following matrices forced to be diagonal even when the original Su is not diagonal (this because DC and PDC/GPDC do not use off-diag terms)
Cd      = diag(diag(Su));   % Cd is useful for calculation of DC
invCd   = inv(Cd);          % invCd is useful for calculation of GPDC
%note: in the definition of the DC here (without inst.eff.) the denominator is not the spectrum because the actual Su is not diagonal

%% Computation of spectral functions

for n = 1:N % at each frequency

    %%% Coefficient matrix in the frequency domain
    As = zeros(M,M); % matrix As(z)=I-sum(A(k))
    for k = 1:p+1
        As = As + A(:,k*M+(1-M:0))*exp(-z*(k-1)*f(n));  %indicization (:,k*M+(1-M:0)) extracts the k-th M*M block from the matrix B (A(1) is in the second block, and so on)
    end

    %%% Transfer matrix (after Eq. 6)
    H(:,:,n)  = inv(As);

    %%% Spectral matrix (Eq. 7)
    S(:,:,n)  = H(:,:,n)*Su*H(:,:,n)'; % ' stands for Hermitian transpose

    %%% Inverse Spectral matrix
    P(:,:,n) = inv(S(:,:,n)); % P(:,:,n) = As'*invSu*As;

    %%% denominators of DC, PDC, GPDC for each m=1,...,num. channels
    for m = 1:M
        tmp1(m)=sqrt((abs(H(m,:,n)).^2) * diag(Cd)); % for the DC: m-th row of H * variance of W (Cd is diag)
        tmp2(m)=sqrt((abs(H(m,:,n)).^2) * ones(M,1)); % for the DTF - don't use covariance information
        tmpp1 = squeeze(As(:,m)); % this takes the m-th column of As...
        tmp3(m) = sqrt(tmpp1'*tmpp1); % for the PDC - don't use covariance information
        tmp4(m) = sqrt(tmpp1'*invCd*tmpp1); % for the GPDC - uses diagonal covariance information
        % for the GPDC - uses diagonal covariance information

    end

    %%% Directed Coherence (Eq. 11)
    DC(:,:,n) = H(:,:,n)*sqrt(Cd) ./ tmp1(:,ones(M,1));
    %nota: tmp1(:,ones(M,1)) crea M colonne tutte uguali a tmp1 - la riga (ossia il den) è la stessa - trova in un colpo solo tutti i denominatori per DC

    %%% Directed Transfer Function (Eq. 11 without sigmas)
    DTF(:,:,n) = H(:,:,n) ./ tmp2(:,ones(M,1));

    %%% Partial Directed Coherence (Eq. 15 without sigmas)
    PDC(:,:,n)  = As./tmp3(ones(1,M),:);
    %nota: tmp3(ones(1,M),:) crea M righe tutte uguali a tmp3 - la colonna (ossia il den) è la stessa - trova in un colpo solo tutti i denominatori per PDC

    %%% Generalized Partial Directed Coherence (Eq. 15)
    GPDC(:,:,n) = (sqrt(invCd)*As) ./ tmp4(ones(1,M),:);

    %%% ISOLATED Coherence
    IPDC(:,:,n) = (sqrt(invCd)*As) ./ sqrt(tmpp1'*invCd*tmpp1);

end

%%% COHERENCE and PARTIAL COHERENCE (Eq. 3)
for m1=1:M
    for m2=1:M
        COH(m1,m2,:) = (S(m1,m2,:)) ./ sqrt(abs(S(m1,m1,:).*S(m2,m2,:)));
        PCOH(m1,m2,:) = (-P(m1,m2,:)) ./ sqrt(abs(P(m1,m1,:).*P(m2,m2,:)));
    end
end

%other definition for partial coherence - Yacoub 1970 - they are really equivalent :)
for n=1:N
    for m1=1:M
        for m2=1:M
            if m1~=m2
                Bb=[P(m1,m1,n) P(m1,m2,n); P(m2,m1,n) P(m2,m2,n);];
                Rb=inv(Bb);
                PCOH2(m1,m2,n) = Rb(1,2) / sqrt( abs(Rb(1,1)) * abs(Rb(2,2)) );
            else
                PCOH2(m1,m2,n)=1;
            end
        end
    end
end
