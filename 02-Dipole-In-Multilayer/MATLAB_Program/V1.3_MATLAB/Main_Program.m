%% *********************************************
tic
showtext=strcat(datestr(now,	'yyyy-mm-dd HH:MM:SS'),': The Program Begin \n');
fprintf(showtext);

%% *********************************************
Main_Basic;
Main_Green;
Main_CalDifference;


%% *********************************************
filename=['./Data/DifferenceReNor_',datestr(now,	'yy_mmdd_HHMM'),'.mat'];
save(num2str(filename));

time=toc;

showtext=strcat(datestr(now,	'yyyy-mm-dd HH:MM:SS'),': Total Simulation Time Is : ',' ',num2str(time),'s \n');
fprintf(showtext);

showtext=strcat(datestr(now,	'yyyy-mm-dd HH:MM:SS'),': Program End\n');
fprintf(showtext);

exit;


