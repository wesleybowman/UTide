clear all

ts=735604;
duration=35;
%duration=360*10;
dt=1/24.0; %in days
time=double(ts:dt:ts+duration);

time_origin=datenum('1899-12-31 12:00:00');% -0.69/24;

time=time;%-time(1)+time_origin;

load('ut_constants.mat','const','shallow');

%period=1/const.freq(48)*3600
%period = 44714.16;
amp =1.0;

phase = 53;

lat=45.5;

period=1./const.freq*3600;
%period(39) = 44714;

%for jj =37:length(period)
for jj =48:48

%time_series=amp*cos((time)/(period(jj)/(24*3600))-phase*pi/180);


time_series=amp*cos((((time-time_origin)*(2*pi/period(jj))*(24*3600))-2*pi*phase/double(360)));

coef = ut_solv(time,time_series ,[] , lat, 'auto','GwchNone','NodsatNone','NoTrend','Rmin',0.95,'OLS','NoDiagn','LinCI');

save ../matlabcoef.mat coef
%[nameu,fu,tidecon,xout]=t_tide(time_series,'starttime',time(1),'latitude',lat);

amp_err(jj)=amp-coef.A(1);

phase_err(jj)=phase-coef.g(1);

ts_recon=ut_reconstr(time, coef);

err(jj)= sqrt(mean((time_series-ts_recon).^2));

ts_fvcom=coef.A(1)*cos(2*pi*((time-mean(time))/(period(jj)/(24*3600))-coef.g(1)/360));

figure
plot(time,time_series)
hold on
plot(time,ts_recon,'k')
plot(time,ts_fvcom,'r')
hold off
datetick

end


% figure
% plot(amp_err)
% hold on
% plot(phase_err,'r')
% plot(err,'k')
% hold off
