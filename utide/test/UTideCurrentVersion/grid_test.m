clear all


filename='/home/wesley/ncfiles/smallcape_force_0001.nc';
ncid=netcdf.open(filename,'NC_NOWRITE');

x = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'x'));
y = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'y'));
lon = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon'));
lat = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat'));
lonc = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lonc'));
latc = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'latc'));
time = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time'));
ua = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ua'),[0 0],[1 length(time)])';
va = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'va'),[0 0],[1 length(time)])';
netcdf.close(ncid);

time=time+datenum('17-Nov-1858');  %On this data, the MJD = 0

order={'M2','S2','N2','K2','K1','O1','P1','Q1'}

coef = ut_solv(time, double(ua), double(va), latc(1), 'auto','NoTrend','Rmin',0.95,'OLS','NoDiagn','LinCI');

save ../speedmatlabcoef.mat coef

[u, v] = ut_reconstr(time, coef);
save ../speedmatlabrecon.mat u v

coef = ut_solv(time, double(ua), [], lat(1), 'auto','NoTrend','Rmin',0.95,'OLS','NoDiagn','LinCI');

save ../matlabcoef.mat coef

[ts_recon, ~]=ut_reconstr(time, coef);
save ../matlabrecon.mat ts_recon
