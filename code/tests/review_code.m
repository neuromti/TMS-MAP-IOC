%% Review Code
mlintrpt('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\+utils','dir')
mlintrpt('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\','dir')
mlintrpt('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\scripts','dir')
%% Check Dependencies
PList = {};
FList = {};
[fList,pList] = matlab.codetools.requiredFilesAndProducts('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\scripts\config_study.m');
PList = [PList,pList.Name]; FList = [FList,fList];

[fList,pList] = matlab.codetools.requiredFilesAndProducts('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\scripts\mapping_study_processing.m');
PList = [PList,pList.Name]; FList = [FList,fList];

[fList,pList] = matlab.codetools.requiredFilesAndProducts('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\scripts\ioc_study_processing.m');
PList = [PList,pList.Name]; FList = [FList,fList];

[fList,pList] = matlab.codetools.requiredFilesAndProducts('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\+utils');
PList = [PList,pList.Name]; FList = [FList,fList];
%%

D1 = dir('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\+utils\*.m');
D2 = dir('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\*.m');
D3 = dir('C:\PROJECTS\Subject Studies\TMS-MAP-IOC\code\scripts\*.m');
D  = cat(1,D1,D2,D3);
PackList = {D.name}';

UFList = unique(FList);
MList = {};
UsedList = {};
% Excluding built-in and fieldtrip functions
for u_idx = 1 : length(UFList),
    a = regexp(UFList(u_idx),('toolbox+.+local'));
    b = regexp(UFList(u_idx),('fieldtrip'));
    if (isempty(a{1}) && isempty(b{1}))
        MList = cat(1,MList,UFList(u_idx));
        [~,fooname,fooext] = fileparts(UFList{u_idx});
        UsedList = cat(1,UsedList,[fooname,fooext]);
    end
end
clc

mlint(MList)
%%
clc
disp(table(unique(PList)','VariableNames',{'Toolbox_Requirements'}))
disp(table(UsedList(~ismember(UsedList,PackList)),'VariableNames',{'Unpackaged'}))
disp(table(PackList(~ismember(PackList,UsedList)),'VariableNames',{'Uncalled'}))
