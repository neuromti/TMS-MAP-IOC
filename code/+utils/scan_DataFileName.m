function T = scan_DataFileName(filename,setup)

    pos.s                     = regexp(filename,'S\w*');
    pos.c                     = regexp(filename,'C\w*');
    pos.d                     = regexp(filename,'.mat');
    T.data_sub              = int32(str2double(filename(pos.s+1:pos.c-1)));
    T.data_cond             = int32(str2double(filename(pos.c+1:pos.d-1)));
    T.subID                 = find(setup.SUB.id==T.data_sub);   
    
end