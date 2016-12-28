% cd('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\')
% clc
% results = runtests('tests\iocdata2statisticsTEST.m');
%%
function tests = iocdata2statisticsTEST

    import matlab.unittest.constraints.IsEqualTo
    
    % Global variables, as local testfunction may have only one input argument
    global TestResults
    
    % Calculate Results
     addpath('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code');
    load('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\tests\testconfigs.mat','DESIGN','STIM','SUBID','SIMPLE_DESIGN');
    DATA = normrnd(0,.5,size(DESIGN,1),1);
    MODULATOR   = (DESIGN(:,1))+(STIM-1)+(SUBID*.1);
    DATA        = DATA+(5*MODULATOR);
      
    [TestResults,ClusterResults] = utils.iocdata2statistics(DATA,DESIGN,SUBID,STIM,true,true);     
    
    % Look for local test functions
    tests  = functiontests(localfunctions);

end

function testAnalytical(testCase)   
    global TestResults
    true_P      = ones(7,4); true_P(:,1)=0;  
    % Checks  whether tests recovered the true effect    
    measuredNorm    = norm(TestResults.Pval-true_P(:,1:3),2); %minimum norm for recovery of factor influence, and 
    randomNorm      = norm(TestResults.Pval-(ones(7,3)),2); %ninimum norm for evenly distributed p-value (i.e. Ho is true)    
    testCase.verifyEqual(measuredNorm<randomNorm,true,'Analytical Analysis did not recover factor influence')
end
% 
function testPermutation(testCase)
    global TestResults
    true_P      = ones(7,4); true_P(:,1)=0;  
    measuredNorm    = norm(TestResults.PermutationPval-true_P,2);  %minimum norm for recovery of factor influence, and 
    randomNorm      = norm(TestResults.PermutationPval-(ones(7,4)),2); %ninimum norm for evenly distributed p-value (i.e. Ho is true)    
    testCase.verifyEqual(measuredNorm<randomNorm,true,'Permutation Analysis did not recover factor influence')
   
end

function testBootStrap(testCase)
    global TestResults
    testCase.verifyEqual(all(TestResults.Power(:,fct)>.9),true,'Bootstrapped power estimate too low')    
end


%     r = corr(SIMPLE_DESIGN(fct,:)',(squeeze(TestResults.MargMeans(:,:)))');  
%     assert(all(r>0.9),['Correlation between Design and Measurements off at Fct',num2str(fct)])
%           
%         
%     r = corr(SIMPLE_DESIGN(fct,:)',reshape(permute(TestResults.CI,[2,1,3]),8,14));  
%     assert(all(r>0.9),['Correlation between Design and CI estimates off at Fct',num2str(fct)])
%   
    
