function sub = intialize_dataset(filename,setup)

t                           = utils.scan_DataFileName(filename,setup);

% Transfer temporary into sub
sub                         = struct();
sub.subID                   = t.subID;

sub.condition               = t.data_cond;
sub.condition_label         = setup.MAP.label.all{sub.condition};
sub.DesignMatrix            = logical([setup.MAP.BI(sub.condition),setup.MAP.LM(sub.condition)]);
sub.sex                     = setup.SUB.sex(t.subID);
sub.age                     = setup.SUB.age(t.subID); 