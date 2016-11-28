function sub = intialize_dataset(filename,setup)

t.s                     = regexp(filename,'S\w*');
t.c                     = regexp(filename,'C\w*');
t.d                     = regexp(filename,'.mat');
t.data_sub              = int32(str2double(filename(t.s+1:t.c-1)));
t.data_cond             = int32(str2double(filename(t.c+1:t.d-1)));
t.subID                 = find(setup.SUB.id==t.data_sub);    

% Transfer temporary into sub
sub                         = struct();
sub.subID                   = t.subID;

sub.condition               = t.data_cond;
sub.condition_label         = setup.MAP.label.all{sub.condition};
sub.DesignMatrix            = logical([setup.MAP.BI(sub.condition),setup.MAP.LM(sub.condition)]);
sub.sex                     = setup.SUB.sex(t.subID);
sub.age                     = setup.SUB.age(t.subID); 