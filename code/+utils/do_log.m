function do_log(logfilename,func_handle)
%takes a logfilename (path\filename) and a anaonymous function handle for a fprintf(logfileid,... to
% and attaches it to the logfile specified in logfilename
logfileid   = fopen(logfilename,'at');
func_handle();
fclose(logfileid);
