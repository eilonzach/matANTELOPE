function epoch_time = serial2epoch(serial_time)
% epoch_time = serial2epoch(serial_time) function to convert serial time,
% in days since 1 Jan 0000 to epochal time in seconds since 1st Jan 1970.
% The former will look like #e5, the latter like #e9

epoch_time = str2epoch(datestr(serial_time,'yyyy-mm-dd HH:MM:SS.FFF'));

end
