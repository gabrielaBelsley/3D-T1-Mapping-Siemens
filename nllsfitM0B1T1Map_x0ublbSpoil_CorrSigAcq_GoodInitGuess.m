function [m0MapNLLS,b1MapNLLS,t1MapNLLS,rsqNLLS,ciT1,NormSSE,deltaT1,Norm2Residuals,ResidualsFit,CorrFactorSpoil] = nllsfitM0B1T1Map_x0ublbSpoil_CorrSigAcq_GoodInitGuess(faArrayRad,faData,sliceVFADixon,TR,maxB1bias,b1CorrectionFactor,lbT1,ubT1,SpoilCorr,varargin)
                                                                                                                            
% nllsfitM0B1T1Map Non Linear-least squares fit on 3D VIBE SPGR data set

%     Inputs:
%       faArrayRad: vector of nominal FAs
%       faData: VFA Data
%       sliceVFADixon: number of SPGR slices to calculate T1
%       TR: Repetition Time in ms
%       maxB1bias: max B1+ factor value 
%       b1CorrectionFactor: B1+ factor 
%       lbT1,ubT1: lower and upper T1 bound used in the incomplete spoiling
%       simulation in ms
%       SpoilCorr: Incomplete Spoiling correction factor

%     Outputs: All outputs at the VFA resolution
%       m0MapNLLS: M0 Map
%       b1MapNLLS: B1+ Map
%       t1MapNLLS: T1 Map
%       rsqMapNLLS: R squared Map
%       ciT1MapNLLS: T1 confidence Interval
%       NormSSE: Normalized Sum Squared Errors: Sum Squared Error for each pixel normalized by the mean signal acquired for the nFA_VFA FAs
%       deltaT1: difference in the upper and lower bound of the T1 Confidence interval
%       Norm2Residuals: Squared norm of the residual
%       ResidualsFit: Residuals
%       CorrFactorSpoil: Incomplete Spoiling correction factor

% Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021

if isempty(sliceVFADixon)
    sliceT1 = 1:size(faData,1);
    nslices = size(faData,1);
    
else
    sliceT1 =sliceVFADixon;
    nslices = length(sliceVFADixon);
end


[nrowsVIBE,ncolsVIBE]= size(faData{sliceVFADixon(1,1),1});

% Intialize variables
DataPerSlicePerPixel = zeros(1,length(faArrayRad));
m0MapNLLS = zeros(nrowsVIBE,ncolsVIBE,nslices);
b1MapNLLS = zeros(nrowsVIBE,ncolsVIBE,nslices);
t1MapNLLS = zeros(nrowsVIBE,ncolsVIBE,max(sliceVFADixon));
rsqNLLS = zeros(nrowsVIBE,ncolsVIBE,nslices);
ciT1 = cell(nrowsVIBE,ncolsVIBE,nslices);
NormSSE = zeros(nrowsVIBE,ncolsVIBE,nslices);
deltaT1 = zeros(nrowsVIBE,ncolsVIBE,nslices);
ResidualsFit = cell(nrowsVIBE,ncolsVIBE,nslices);
Norm2Residuals = zeros(nrowsVIBE,ncolsVIBE,nslices);
CorrFactorSpoil = zeros(nrowsVIBE,ncolsVIBE,nslices,length(faArrayRad));

% SPGR function used to calculate M0
funSPGR = @(TR,alpha,T1)((sin(alpha).*(1-exp(-TR/T1)))./(1-cos(alpha).*exp(-TR/T1)));

tol = 1; % stop iteration when change between T1 prefit and after fit changes by less than 1 ms


%---22/01/2020: Trust-region-reflective and levenberg-marquardt gives
% identical T1 values (difference on the order of e-4 or less). Because we
% know the boundaries of the liver T1 we want to map (600-1000ms) we use
% trust region reflective as levenberg-marquardt does not allow for bounds
% and thus takes longer/more iterations.
%Option1: Trust-Region-Reflective
options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','off','OptimalityTolerance',1.0000e-10);
%Option 2: Levenberg-Marquardt
%options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','Display','off');

%----Trust region reflective NLLS algorithm------
%https://www.mathworks.com/help/optim/ug/lsqnonlin.html#buuhch7-options
%https://www.mathworks.com/help/optim/ug/least-squares-model-fitting-algorithms.html#f204
%https://www.mathworks.com/help/optim/nonlinear-programming.html?s_tid=CRUX_lftnav
%https://www.mathworks.com/help/optim/ug/optimization-toolbox-tutorial.html
%https://www.mathworks.com/help/optim/ug/problem-based-workflow.html
%https://www.mathworks.com/help/optim/ug/example-nonlinear-constrained-minimization.html#brg0p3g-9
%----Summary Optimisation algorithms------
%https://www.mathworks.com/help/optim/ug/optimization-decision-table.html
%----Analysing the Output of nlls: number iterations, Stopping Criteria------
% https://www.mathworks.com/help/optim/ug/when-the-solver-might-have-succeeded.html#br544l9-1
%https://www.mathworks.com/help/optim/ug/iterative-display.html
% https://www.mathworks.com/help/optim/ug/first-order-optimality-measure.html
% https://www.mathworks.com/help/optim/ug/tolerances-and-stopping-criteria.html
% https://www.mathworks.com/help/optim/ug/optimizing-a-simulation-or-ordinary-differential-equation.html

for islice = sliceT1
    fprintf('Calculating T1 Slice %i \n', islice)

   
    for irow = 1:1:nrowsVIBE
       
        for icol = 1:1:ncolsVIBE
            

 
            % Only calculate T1 if B1+ is within the range we optimised the B1+ precision
            if (b1CorrectionFactor(irow,icol,islice) > 0.1) && (b1CorrectionFactor(irow,icol,islice) < maxB1bias)
                
                % Read the data for each pixel
                for ifa = 1:length(faArrayRad)
                    DataPerSlice = double(faData{islice,ifa});
                    DataPerSlicePerPixel(1,ifa) = DataPerSlice(irow,icol);
                    
                end
                

                % True FA: Apply the B1+ Correction to the nominal FA
                    faTrueRad = faArrayRad.*b1CorrectionFactor(irow,icol,islice);
                    faTrueDeg = (b1CorrectionFactor(irow,icol,islice).*rad2deg(faArrayRad));
                    
                    [T1_InitialEst] = nlls_T1InitialEstimate(DataPerSlicePerPixel,faTrueRad,TR,lbT1,ubT1);
 
                % if the Signal for any FA is zero or very low(50/2^12)
                % it's not worth fitting the data and we place t1MapNLLS to zero
                if any(DataPerSlicePerPixel(1,:)== 0) || any(DataPerSlicePerPixel(1,:)<20) 

                           t1MapNLLS(irow,icol,islice) = 0;
                    
                    
                else
                    
                    T1_PreFit = 10000; %large number such that the first while argument is always true
                    T1_AfterFit = T1_InitialEst; % Initial T1 estimate 
                    iter = 0;
                    
                    %if the difference between the previous T1 estimate
                    %(T1_PreFit) and the one just estimated T1_AfterFit is
                    %more than 1ms (tol) then we do another fit
                    while abs(T1_AfterFit - T1_PreFit)>tol && iter<5 
                        T1_PreFit = T1_AfterFit;
                        iter = iter+1;
                       
                        %-------Spoiling correction--------------------------------------------------------------------
                        %Default T1 for SpoilCorr is T1_PreFit
                        T1_SpoilCorr = T1_PreFit;
                        % Verify it is not outside the T1 bounds for which the
                        % Spoil Correction was calculated. 
                        lbT1Spoil = lbT1;
                        ubT1Spoil = ubT1;
                         
                        if T1_PreFit <  lbT1Spoil
                            T1_SpoilCorr = lbT1Spoil;
                        elseif T1_PreFit >  ubT1Spoil
                            T1_SpoilCorr = ubT1Spoil;
                        end
                            
                        % Spoiling correction uses the trueFA
                        if isempty(SpoilCorr)
                            DataPerSlicePerPixel_SpoilCorr= DataPerSlicePerPixel;
                        else
                            CorrFactorSpoil(irow,icol,islice,:) = SpoilCorr(repmat(T1_SpoilCorr,[1,length(faTrueDeg)]),faTrueDeg);
                            DataPerSlicePerPixel_SpoilCorr = DataPerSlicePerPixel.*(squeeze(CorrFactorSpoil(irow,icol,islice,:)).');
                        end
                        %---------------------------------------------------------------------------------------------
                        
                        if isempty(varargin) %estimate M0 and T1
                            %use the smallest angle faArrayRad(1,1) and
                            %signal DataPerSlicePerPixel_SpoilCorr(1,1) for an initial estimate of 
                            %M0 as this will be less affected by spoiling
                            M0InitialGuess = DataPerSlicePerPixel_SpoilCorr(1,1)./funSPGR(TR,faTrueRad(1,1),T1_PreFit);
                            x0 = [M0InitialGuess,T1_PreFit];%intial guess of M0, T1
                                                                       
                            % M0 lower bound: subtract 3000 from M0 initial guess, if negative set it to 0
                                M0lb = max(0,M0InitialGuess-3000); 
          
                            
                        else % fix M0 and estimate only T1
                            M0Map = varargin{1,1}; %M0 passed as input
                            x0 = T1_PreFit;%intial guess of T1
                        end
                        
                        
                        if length(x0) == 3 % solve for M0, T1, B1
                            lb = [M0lb,lbT1,0.1];                 % lower bound M0, T1, B1
                            ub = [M0InitialGuess+20000,ubT1,maxB1bias]; % upper bound M0, T1, B1
                            
                            fun = @(x)(x(1)*(sin(faArrayRad.*x(3)).*(1-exp(-TR/x(2))))./(1-cos(faArrayRad.*x(3)).*exp(-TR/x(2))))-DataPerSlicePerPixel_SpoilCorr;
                            % x is a vector with three elements: x(1)=M0; x(2) = T1; x(3) = B1Corr
                            
                        elseif length(x0) == 2 % solve for M0, T1
                            lb = [M0lb,lbT1];                 % lower bound M0, T1
                            ub = [M0InitialGuess+20000,ubT1]; % upper bound M0, T1
                            
                            % note: use the fa in rad after applying B1+ map (faTrueRad) and
                            % the spoil corrected signal
                            fun = @(z)(z(1)*(sin(faTrueRad).*(1-exp(-TR/z(2))))./(1-cos(faTrueRad).*exp(-TR/z(2))))-DataPerSlicePerPixel_SpoilCorr;
                            

                        elseif length(x0) == 1 % solve for T1 only, fix M0
                            
                            lb = lbT1; % lower bound T1
                            ub = ubT1; % upper bound T1
  
                            fun = @(x)(M0Map(irow,icol,islice)*(sin(faArrayRad.*b1CorrectionFactor(irow,icol,islice) ).*(1-exp(-TR/x(1))))./(1-cos(faArrayRad.*b1CorrectionFactor(irow,icol,islice)).*exp(-TR/x(1))))-DataPerSlicePerPixel_SpoilCorr;
                            
                        end
                        
                        % Note on algorithm options: https://uk.mathworks.com/help/optim/ug/choosing-the-algorithm.html#bsbwx1l
                        %try 'trust-region-reflective' first. If your problem has bounds, you must use 'trust-region-reflective'.
                        % If your problem has no bounds and is underdetermined (fewer equations than dimensions), use 'levenberg-marquardt'
                        
                        
                        % Solve the non-linear least squares fit using
                        % trust-region reflective, initial guess and upper and lower bounds
                        [x,resnorm,residuals,~,~,~,J] = lsqnonlin(fun,x0,lb,ub,options);
                        %[x,resnorm,residual,exitflag,output,lambda,jacobian]
                        

                        if length(x0) == 3
                            m0MapNLLS(irow,icol,islice) = x(1);
                            t1MapNLLS(irow,icol,islice) = x(2);
                            b1MapNLLS(irow,icol,islice) = x(3);
                        elseif length(x0) == 2
                            m0MapNLLS(irow,icol,islice) = x(1);
                            
                            t1MapNLLS(irow,icol,islice) = x(2);
                            ResidualsFit{irow,icol,islice} = residuals;
                            Norm2Residuals(irow,icol,islice) = resnorm;

                            
                        elseif length(x0) == 1
                            t1MapNLLS(irow,icol,islice) = x(1);
                        end
                        

                        % Note: to estimate the CI one needs at least one more
                        % data point than the number of variables to estimate
                      ciT1{irow,icol,islice} = nlparci(x,residuals,'jacobian',J);
                      
                      if ~isempty(ciT1{irow,icol,islice})
                          deltaT1(irow,icol,islice) = ciT1{irow,icol,islice}(2,2)- ciT1{irow,icol,islice}(2,1);
                      end
                        
                        y_signal = DataPerSlicePerPixel_SpoilCorr;
                       rsqNLLS(irow,icol,islice) = 1 - ((sum(fun(x)).^2)/sum((y_signal - mean(y_signal)).^2));
                        % Average Fractional Error Per Angle: Sum Squared Error for each pixel normalized by
                        % the mean signal acquired for the 5FA
                       NormSSE(irow,icol,islice) = sqrt(resnorm)/(mean(y_signal).*length(faArrayRad));
                        %normalized by the mean(y_signal): so that we can compare with Trio and Prisma
                        %normalized by the length(faArrayRad): so that we can compare between
                        %different number of FAs

                        clear x;
                   
                        T1_AfterFit = t1MapNLLS(irow,icol,islice);
                        
                       %  fprintf('T1 = %10.2f Iter = %i\n', T1_AfterFit, iter)
                    
                    
                    end
                    
                end
                
            end

                
        end
    end
end



