function serial_time = epoch2serial(epoch_time)
% serial_time = epoch2serial(epoch_time)
% function to convert epochal time in seconds since 1st Jan 1970 to serial
% time, in days since 1 Jan 0000. The former will look like #e9, the latter
% like #e5

serial_time = datenum(epoch2str(epoch_time,'%Y-%m-%d %H:%M:%S'));

end
