%% load all the data
%[timearray, dataarray, dataheader] = abiss_parse_log_file(fullpath)
datestamp='20180915';
timestamp='103549';
timestamp='133434';
timestamp='174941';
dirpath = [ '/Users/hollieemery/Dropbox/Harvard/ABISS/RV Falkor (September 2018; McArthur Ridge)/2018 Deployment Data'];
[to1a, o1a, o1h] = abiss_parse_log_file( [dirpath filesep 'abiss_1' datestamp '_' timestamp '.dat' ] );
[to2a, o2a, o2h] = abiss_parse_log_file( [dirpath filesep 'abiss_2' datestamp '_' timestamp '.dat' ] );
[tva, va, vh] = abiss_parse_log_file( [dirpath filesep 'abiss_V' datestamp '_' timestamp '.dat' ] );
[tra, ra, rh] = abiss_parse_log_file( [dirpath filesep 'abiss_R' datestamp '_' timestamp '.dat' ] );

%% make timetables for export 
o2time_1 = seconds((to1a-min(to1a)));
o2time_2 = seconds((to2a-min(to2a)));
ra_time = seconds((tra-min(tra)));
va_time = seconds((tva-min(tva)));

o2_1 = table2timetable(cell2table(o1a,'VariableNames',o1h),'RowTimes',o2time_1);
o2_2 = table2timetable(cell2table(o2a,'VariableNames',o2h),'RowTimes',o2time_2);
ra = table2timetable(table(ra,'VariableNames',rh),'RowTimes',ra_time);

%% plot battery log
plot(tva, va)

%% export to csv
writetimetable(o2_1,'o2_1.csv')
writetimetable(o2_2,'o2_2.csv')
writetimetable(ra,'ra.csv')
