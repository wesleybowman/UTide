function coef = ut_solv(tin,uin,vin,lat,cnstit,varargin)
% UT_SOLV() Execute harmonic tidal analysis.
%   Companion function UT_RECONSTR() reconstructs superposed harmonics
%       for hindcast and/or forecast/prediction.
%
% OVERVIEW 
% 
% Syntax for two-dimensional raw input, such as velocities:
% 	coef = UT_SOLV ( t_raw, u_raw, v_raw, lat, cnstit , {options} ); 
%   [ u_fit, v_fit ] = UR_RECONSTR ( t_fit, coef , {options} );
% 
% Syntax for one-dimensional raw input, such as sea level:
% 	coef = UT_SOLV ( t_raw, sl_raw, [], lat, cnstit , {options} ); 
%   [ sl_fit, ~ ] = UT_RECONSTR ( t_fit, coef , {options} ); 
% 
% Input/output parameters (more detailed explanations below): 
%   * t_raw = times of raw inputs
%   * u_raw & v_raw, or sl_raw = raw input values
%   * lat = latitude (required for use in default nod./sat. corrections)
%   * cnstit = specification of constituents to include in model
%   * coef = results structure from UT_SOLV(), used by UT_RECONSTR()
%   * t_fit = arbitrary times for reconstructed superposed harmonics
%   * u_fit & v_fit, or sl_fit = reconstructed superposed harmonics
%   * {options} explained below
%
% Analysis of groups of records
%   The descriptions that follow next for INPUTS, DEFAULTS, OPTIONS, and 
%   OUTPUT are for treatment of a single record. Following that, 
%   explanations are given for modifications that enable treating a group 
%   of records with a single execution.
%
% INPUTS
% 
% t_raw
%
%   * Column vector of datenum values for the sampled times, UTC/GMT.
%   * May have NaNs; they are removed during analysis, along with the 
%       corresponding (NaN or non-NaN) u_raw/v_raw or sl_raw values.
%   * May be uniformly distributed ("equispaced") or irregular.
%   * Considered equispaced if (after NaNs are removed) 
%           var(unique(diff(t_raw)))<eps 
%       is true, in which case FFT methods [ pwelch(), cpsd() ] are used 
%       (after filling any NaNs in u_raw/v_raw or sl_raw by interpolation,
%       as needed) for the periodogram of the residual in the colored 
%       confidence interval calculation.
%   * Considered irregularly distributed if (after NaNs are removed)
%           var(unique(diff(t_raw)))>=eps 
%       is true, in which case the Lomb-Scargle periodogram of the 
%       residual is used in the colored confidence interval calculation. 
%
% u_raw and v_raw, or sl_raw
%
%   * Real-valued column vectors of the same size as t_raw.
%   * May have NaNs; they are removed during analysis.
%
% lat
%
%   * Latitude (in decimal degrees, positive North).
%   * Ignored if nodal/satellite corrections omitted ('NodsatNone').
% 
% cnstit 
%
%   * Specification of constituents included in model.
%   * One of the following:
%       The string ‘auto’ (not case-sensitive), to implement the F77
%           automated decision tree (default Rmin = 1). If inference/
%           reference constituents are also specified they are included 
%           regardless of the decision tree results.
%       A cell array of constituent names (4-character strings, not 
%           case-sensitive). Constituents available (including
%           shallow-water constituents) are those in the const.name 
%           variable in the ut_constants.mat file. 
%
% {options}
%
%   * See section below.
%
% DEFAULTS
%
% Defaults (when no {options} parameters are passed in):
%   * The linear (secular, non-tidal) trend is included in the model.
%   * No pre-filtering correction is made.
%   * Nodal/satellite corrections with exact times are implemented.
%   * Greenwich phase lags use the astronomical argument with exact times.
%   * No constituents are inferred.
%   * If cnstit is ‘auto’, the automated decision tree constituent 
%       selection method is applied with Rmin = 1.
%   * The solution method is iteratively re-weighted least squares (IRLS)
%       with Cauchy weight function and tuning parameter 2.385, the Matlab
%       default value for Cauchy (TunRdn = 1). See Matlab robustfit() help. 
%   * Confidence intervals use Monte Carlo, 200 realizations (Nrlzn=200).
%   * Confidence intervals use (colored) spectra of actual residuals.
%   * For irregular input times, spectral calculations use Lomb-Scargle 
%       periodograms with frequency oversampling factor 1 (LSFrqOSmp=1).
%   * The constituent selection diagnostics table is generated, and uses
%       minimum signal-to-noise (SNR) of MinSNR = 2.
%   * The order of values in the fields of coef and coef.diagn that have 
%       size nallc x 1 is by decreasing percent energy (PE) of the
%       corresponding constituents.
%   * During runtime the following are displayed: (1) results 
%       (coefficients with confidence intervals), (2) run description 
%       meta information, and (3) constituent selection diagnostics table. 
% 
% OPTIONS
%
% Option flags are not case-sensitive but cannot be abbreviated. The order
%   of the option flags is not important but they must be passed in 
%   after all other arguments. As shown below, some options require an 
%   associated parameter. See report for more complete explanations.
%
%   ‘NoTrend’
%
%       * Omit the linear/secular trend term from the model. 
%
%   ‘PreFilt’, PreFilt 
%
%       * Implement the correction for pre-filtering applied to raw inputs.
%           PreFilt.P = np x 1  vector of P values.
%           PreFilt.frq = np x 1 vector of frequencies (cycles per hour). 
%           PreFilt.rng = two-value vector with min/max allowable P. 
%
%   One of the following:
%
%       ‘NodsatLinT’ 
%           * Use linearized times in nodal/satellite corrections.
%
%       'NodsatNone'
%           * Omit nodal/satellite corrections.
%
%   One of the following:
%
%       'GwchLinT'
%           * Use linearized times in astronomical argument when computing
%               Greenwich phase lags.
%
%       ‘GwchNone’
%           * Report raw (not Greenwich-referenced) phase lags.
%
%   ‘Infer’,Infer
%
%       * Include nI inference constituents in model.
%           Infer.infnam = cell array (nI x 1) of 4-character names 
%               of inferred constituents.
%           Infer.refnam = same, for reference constituents; need not 
%               all be unique unless using 'InferAprx'.
%           Infer.amprat = real-valued amplitude ratios (unitless)
%           Infer.phsoff = real-valued phase offsets (degrees)
%               For 2D case, Infer.amprat and Infer.phsoff are each 
%                   2nI x 1 with the + (-) values in the first (last) nI 
%                   elements; in the 1D case they are nI x 1.
%
%   ‘InferAprx’
%
%       * Implement approximate method for inference calculation. Ignored 
%           unless ‘Infer’ also selected. 
%
%   ‘Rmin’, Rmin
%
%       * Specify Rmin value (minimum conventional Rayleigh criterion) to 
%           be used with the automated constituent selection. Default is 1. 
%           Ignored if cnstit is not ‘auto’. 
%
%   One of the following:
%
%       ‘OLS’ 
%           * Use ordinary least squares solution method.
%
%       ‘Andrews’,‘Bisquare’,‘Fair’,‘Huber’,‘Logistic’,‘Talwar’, OR ‘Welsch’ 
%           * Weight function for IRLS method. Default is 'Cauchy'. See
%               Matlab documentation for robustfit().
%
%   ‘TunRdn’, TunRdn 
%
%       * Reduce IRLS tuning parameter relative to Matlab default, for 
%           the given weight function, by the tuning parameter reduction 
%           factor TunRdn (tuning parameter used is default divided by 
%           TunRdn). Default is TunRdn = 1. Ignored if using ‘OLS’.
%
%   ‘LinCI’ 
%
%       * Use linearized method to determine confidence intervals.
%
%   ‘White’ 
%
%       * Use white noise floor assumption for confidence intervals.
%
%   ‘Nrlzn’, Nrlzn
%
%       * Use Nrlzn realizations in Monte Carlo confidence interval 
%           calculation. Default Nrlzn=200. Ignored if ‘LinCI’.
%
%   ‘LSFrqOSmp’, LSFrqOSmp
%
%       * Frequency oversampling factor in Lomb-Scargle. Default 1.
%
%   ‘DiagnMinSNR’, DiagnMinSNR
%
%       * Specify minimum SNR for constituents included in the 
%           reconstructed fit used for TVsnrc and PTVsnrc diagnostics, and 
%           in the diagnostic plots (if generated). Default value is 2.
%
%   One of the following:
%
%       ‘NoDiagn’
%           * Skip diagnostics computations and diagnostic figures.
%
%       'DiagnPlots'
%           * Generate plots in addition to calculating diagnostics;
%               ignored if 'NoDiagn' is specified.
%
%   ‘OrderCnstit’, CnstitSeq
%
%       * Override default PE-ranked ordering of constituent-based 
%           parameters in output structure coef.
%               CnstitSeq is one of the following:
%                   * ‘snr’, to order by decreasing SNR, OR
%                   * ‘frq’, to order by increasing frequency, OR
%                   * a cell array of 4-character strings that, if it 
%                       differs from cnstit, only differs by the order of 
%                       its rows (allowed only for cnstit not ‘auto’).
%   
%   'RunTimeDisp', RunTimeDisp
%
%       * Default runtime display (RunTimeDisp='yyy') includes display of
%           (1) results coefficients with confidence intervals, (2) run
%           description meta information, and (3) the constituent selection
%           diagnostics table. To suppress all runtime display information,
%           use RunTimeDisp = 'nnn'; to suppress only the first and third
%           elements, use 'nyn'; etc. If 'NoDiagn' is selected, the third
%           character in RunTimeDisp will be ignored (no table exists nor
%           will be displayed).
%
% OUTPUT
%
% coef
%
%   * A structure with all analysis results in various scalar, vector, and 
%       character array fields.
%   * Vector fields of coef are n_allc x 1, where n_allc is the total 
%       number of constituents (non-reference, reference, and inferred) 
%       included in the model, and by default the order of the elements is 
%       by decreasing percent energy (PE); this ordering may be changed, 
%       using 'OrderCnstit', to be decreasing SNR, increasing frequency, 
%       or a user-specified order.
% 
% coef.name
%
%   * Cell array (n_allc x 1) of 4-character names of constituents.
%
% coef.Lsmaj, coef.Lsmin, coef.theta, coef.g (2D case) 
% coef.A, coef.g (1D case) 
%
%   * Real-valued column vectors (n_allc x 1) of current ellipse parameters
%       or amplitude/phase. Lsmaj, Lsmin, A are in raw input units, and
%       theta (current ellipse orientation angle) and g (phase lag) are 
%       in degrees.
%
% coef.Lsmaj_ci, coef.Lsmin_ci, coef.theta_ci, coef.g_ci (2D case) 
% coef.A_ci, coef.g_ci (1D case) 
%
%   * Estimated 95% confidence intervals.
%
% coef.umean, coef.vmean (2D case) 
% coef.mean (1D case)
%
%   * Scalar(s) holding real-valued mean(s) (raw input units) determined 
%       by harmonic analysis.
%
% coef.uslope, coef.vslope (2D case) 
% coef.slope (1D case) 
%
%   * Scalar(s) holding real-valued slope(s) (raw input units per day) of 
%       the linear/secular trend determined by harmonic analysis. Omitted 
%       if trend is not included in model ('NoTrend' option).
%
% coef.results
%
%   * Character array summarizing all above-listed elements of coef.
% 
% coef.aux.rundescr
%   
%   * Character array with a descriptive explanation of the run.
%
% coef.aux.opt
%
%   * Fields containing all option settings and input parameters.
%
% coef.aux.frq
%
%   * Real-valued column vector (n_allc x 1) of frequencies of included 
%       constituents (cycles per hour).
%
% coef.aux.lind
%
%   * Column vector (n_allc x 1) of list indices, in ut_constants.mat, of 
%       included constituents. Used by UT_RECONSTR().
%
% coef.aux.lat
%
%   * Latitude provided as input. Used by UT_RECONSTR().
%
% coef.aux.reftime
%
%   * Reference time (mean of first & last input times), datenum UTC/GMT.
%
% coef.diagn (unless 'NoDiagn' used)
%
%   * Series of fields holding results of diagnostic calculations. Vectors
%       are n_allc x 1 with elements ordered by decreasing PE (regardless 
%       of 'OrderCnstit').
%   * The cell string array coef.diagn.table includes all coef.diagn fields
%       presented in a concise manner for visual inspection within Matlab. 
%   * The report provides complete explanations.
%
% GROUPS OF RECORDS
%
% The group of n_s time sequences to be treated is indexed as an
%   n_d - dimensional array of size n1 x n2 x n3 ... n-n_d, where each n 
%   value gives the size of that dimension of the array and none of the n
%   values may be 1.
%
% Inputs to UT_SOLV() are as in the single-record case except that:
%   
%   * t_raw can be either (1) a single n_t x 1 vector of times that applies 
%       to all members of the group, OR (2) an n_t x n1 x n2 x n3 ... x 
%       n-n_d array of times, when more than one record in the group has a 
%       different set of times; in the latter case the number of times 
%       must be the same (n_t) for each record, which can be accomplished
%       by NaN-padding at the start/end of each shorter record if needed.
%
%   * u_raw and v_raw, or sl_raw, are each n_t x n1 x n2 x n3 ... x 
%       n-n_d arrays.
%
%   * lat is either a scalar (same for all sequences) or an n1 x n2 x n3
%       ... n-n_d array (different value for each sequence).
%
%   * cnstit cannot be 'auto' and must be a specific list of constituents 
%       to include (necessary for convenient grouping of results fields).
%
%   * For inference, (1) the same inference and reference constituents 
%       are used for each record (a single set of Infer.infnam and 
%       Infer.refnam is allowed), and (2) Infer.amprat and Infer.phsoff are 
%       EITHER: 2*nI x 1 (two-dimensional case) or nI x 1 (one-dimensional 
%       case), to apply the same values to every record; OR: 2*nI x n1 x 
%       n2 x n3 ... n-n_d (two-dimensional case) or nI x n1 x n2 x n3 ...
%       x n-n_d (one-dimensional case), to apply different values to 
%       different records.
%
% Options are as in the single-record case described above, and applied 
%   to each member of the group, with the following exceptions: diagnostic 
%   plots are not allowed (no 'DiagnPlots' option); ordering of 
%   constituent-indexed output fields is always by increasing frequency 
%   ('OrderCnstit' is fixed as 'frq'); and there is no runtime display 
%   ('RunTimeDisp' is fixed as 'nnn').
%
% Output is as in the single-record case except for the following changes:
%
%   * Fields of coef and its subfields that are of size n_allc x 1 in the 
%       single-record case have size n_allc x n1 x n2 x n3 ... x n-n_d, 
%       with ordering of the first dimension by increasing frequency.
%
%   * Fields of size 1 x 1 in the single-record case have size n1 x n2 x
%       n3 ... x n-n_d.
%   
%   * Exceptions to the above rules, to avoid redundancies, are:
%       (1) coef.name, coef.aux.frq, and coef.aux.lind are each n_allc x 1.
%       (2) the only fields of coef.aux.opt that potentially differ from 
%           the single-record case are equi, infer.amprat and 
%           infer.phsoff; the latter two are sized as their inputs.
%       (3) if the input lat is a scalar then coef.aux.lat is a scalar.
%
% For more information see:
%   Codiga, D.L., 2011. Unified Tidal Analysis and Prediction Using the 
%       UTide Matlab Functions. Technical Report 2011-01. Graduate School 
%       of Oceanography, University of Rhode Island, Narragansett, RI. 
%       59pp. ftp://www.po.gso.uri.edu/pub/downloads/codiga/pubs/
%       2011Codiga-UTide-Report.pdf
%
% UTide v1p0 9/2011 d.codiga@gso.uri.edu

if isequal(sum(size(uin)>2),1)  % single record
    coef = ut_solv1(tin,uin,vin,lat,cnstit,varargin{:});    
else % group of records
    % n_t = number of times in each sequence in the group
    n_t = size(uin,1); 
    % alln = [n_1 n_2 n_3 ... n_n_d]
    alln = size(uin);
    alln = alln(2:end);
    % n_d = number of dimensions of group
    n_d = length(alln);
    % n_s = number of time sequences in group
    n_s = prod(alln); 
    % size of group
    if isequal(n_d,1)
        gsz = [1 alln];
    else
        gsz = alln;
    end
    % check and condition inputs
    if (~isequal(size(tin,1),size(uin,1)) && ...
            ~isequal(size(tin),size(uin))) || (~isempty(vin) && ...
            ~isequal(size(vin),size(uin)))
        error(['ut_solv(group): sizes of t & u, or of their first '...
            'dimensions, must be equal; u & v sizes must also be equal.']);
    end
    if ~isequal(size(lat),[1 1]) && ~isequal(size(lat),gsz)
        error('ut_solv(group): lat has neither acceptable size.');
    end
    if isequal(lower(cnstit),'auto')
        error(['ut_solv(group): ''cnstit'' cannot be ''auto'' when '... 
            'analyzing groups of records.']);
    end
    ind = [];
    for i = 1:length(varargin)
        if ischar(varargin{i})
            if isequal(lower(varargin{i}),'ordercnstit');
                ind = i;
            end
        end
    end
    if ~isempty(ind)
        if ~isequal(varargin{ind+1},'frq')
            error(['ut_solv(group): if ''OrderCnstit'' is passed in '...
                'it must be ''frq''. ']);
        end
        varargin{ind:ind+1} = [];
    end
    ind = [];
    for i = 1:length(varargin)
        if ischar(varargin{i})
            if isequal(lower(varargin{i}),'runtimedisp');
                ind = i;
            end
        end
    end
    if ~isempty(ind)
        if ~isequal(varargin{ind+1},'nnn')
            error(['ut_solv(group): if ''RunTimeDisp'' is passed in '...
                'it must be ''nnn''. ']);
        end
        varargin{ind:ind+1} = [];
    end
    for i = 1:length(varargin)
        if ischar(varargin{i})
            if isequal(lower(varargin{i}),'diagnplots');
                 error('ut_solv(group): ''DiagnPlots'' is not allowed.');
            end
        end
    end
    infind = [];
    for i = 1:length(varargin)
        if ischar(varargin{i})
            if isequal(lower(varargin{i}),'infer');
                 infind = i;
            end
        end
    end
    if ~isempty(infind)
        infer0 = varargin{infind+1};
        ninf = length(infer0.infnam);
        if ~isempty(vin)
            if ( ~isequal(size(infer0.amprat),[2*ninf 1]) && ...
                    ~isequal(size(infer0.amprat),[2*ninf alln]) ) || ...
                    ( ~isequal(size(infer0.phsoff),[2*ninf 1]) && ...
                    ~isequal(size(infer0.phsoff),[2*ninf alln]) )
                error(['ut_solv(group): infer.amprat and/or '...
                    'infer.phsoff have neither of the acceptable sizes.']);
            end
        else
            if ( ~isequal(size(infer0.amprat),[ninf 1]) && ...
                    ~isequal(size(infer0.amprat),[ninf alln]) ) || ...
                    ( ~isequal(size(infer0.phsoff),[ninf 1]) && ...
                    ~isequal(size(infer0.phsoff),[ninf alln]) )
                error(['ut_solv(group): infer.amprat and/or '...
                    'infer.phsoff have neither of the acceptable sizes.']);
            end
        end
    end
    nttst = 0;
    for i = 1:length(varargin)
        if ischar(varargin{i})
            if isequal(lower(varargin{i}),'notrend')
                nttst = 1;
            end
        end
    end
    ndtst = 0;
    for i = 1:length(varargin)
        if ischar(varargin{i})
            if isequal(lower(varargin{i}),'nodiagn')
                ndtst = 1;
            end
        end
    end
    % reshape group to vectors (from mult-dim arrays) if needed
    if n_d > 1
        if isequal(size(tin),size(uin))
            tin = reshape(tin,n_t,n_s);
        end
        uin = reshape(uin,n_t,n_s);
        if ~isempty(vin)
            vin = reshape(vin,n_t,n_s);
        end
        if ~isempty(infind)
            if ~isempty(vin)
                if ~isequal(size(infer0.amprat),[2*ninf 1])
                    infer0.amprat = reshape(infer0.amprat,2*ninf,n_s);
                    infer0.phsoff = reshape(infer0.phsoff,2*ninf,n_s);
                end
            else
                if ~isequal(size(infer0.amprat),[ninf 1])
                    infer0.amprat = reshape(infer0.amprat,ninf,n_s);
                    infer0.phsoff = reshape(infer0.phsoff,ninf,n_s);
                end
            end
        end
        if ~isequal(size(lat),[1 1])
            lat = reshape(lat,1,n_s);
        end
    end
    % initialize storage (as vectors)
    if ~isempty(infind)
        nallc = length(unique([cnstit; infer0.infnam; infer0.refnam;]));
    else
        nallc = length(cnstit);
    end
    coef.g = nan*ones(nallc,n_s);
    coef.g_ci = coef.g;
    if ~isempty(vin)
        coef.Lsmaj = coef.g;
        coef.Lsmaj_ci = coef.g;
        coef.Lsmin = coef.g;
        coef.Lsmin_ci = coef.g;
        coef.theta = coef.g;
        coef.theta_ci = coef.g;
    else
        coef.A = coef.g;
        coef.A_ci = coef.g;
    end
    if ~isempty(vin)
        coef.umean = nan*ones(1,n_s);
        coef.vmean = coef.umean;
        if ~nttst
            coef.uslope = coef.umean;
            coef.vslope = coef.vmean;
        end
    else
        coef.mean = nan*ones(1,n_s);
        if ~nttst
            coef.slope = coef.mean;
        end
    end
    coef.results = cell(1,n_s);
    coef.aux.rundescr = cell(1,n_s);
    coef.aux.opt.equi = nan*ones(1,n_s);
    coef.aux.reftime = nan*ones(1,n_s);
    if ~ndtst
       coef.diagn.PE = coef.g;
       coef.diagn.SNR = coef.g;
       coef.diagn.TVraw = nan*ones(1,n_s);
       coef.diagn.TVallc = nan*ones(1,n_s);
       coef.diagn.TVsnrc = nan*ones(1,n_s);
       coef.diagn.PTVallc = nan*ones(1,n_s);
       coef.diagn.PTVsnrc = nan*ones(1,n_s);
       coef.diagn.SNRallc = nan*ones(1,n_s);
       coef.diagn.K = nan*ones(1,n_s);
       coef.diagn.lo.CorMx = coef.g;
       coef.diagn.hi = coef.diagn.lo;
       coef.diagn.table = cell(1,n_s); 
    end
    % main loop, over sequences in group
    for is = 1:n_s
        % prep inputs to ut_solv1
        uin1 = uin(:,is);
        if ~isempty(vin)
            vin1 = vin(:,is);
        else
            vin1 = [];            
        end
        if isequal(size(uin),size(tin))
            tin1 = tin(:,is);
        else
            tin1 = tin;
        end
        if ~isequal(size(lat),[1 1])
            lat1 = lat(is);
        else
            lat1 = lat;
        end
        if ~isempty(infind)
            infer1 = infer0;
            if ~isempty(vin) 
                if ~isequal(size(infer0.amprat),[2*ninf 1])
                    infer1.amprat = infer0.amprat(:,is);
                    infer1.phsoff = infer0.phsoff(:,is);
                end
            else
                if ~isequal(size(infer0.amprat),[ninf 1])
                    infer1.amprat = infer0.amprat(:,is);
                    infer1.phsoff = infer0.phsoff(:,is);
                end
            end
            varargin{infind+1} = infer1;
        end
        % execute ut_solv1
        fprintf(sprintf('[%d/%d]',is,n_s));
        try
            coef1 = ut_solv1(tin1,uin1,vin1,lat1,cnstit,...
                'OrderCnstit','frq','RunTimeDisp','nnn',varargin{:});
        catch err
            disp(err.message);
        end
        % store results from 1st run that apply to all runs
        if is == 1
            coef.name = coef1.name;
            coef.aux.frq = coef1.aux.frq;
            coef.aux.lind = coef1.aux.lind;
            coef.aux.opt.twodim = coef1.aux.opt.twodim;
            coef.aux.opt.notrend = coef1.aux.opt.notrend;
            coef.aux.opt.prefilt = coef1.aux.opt.prefilt;
            coef.aux.opt.nodsatlint = coef1.aux.opt.nodsatlint;
            coef.aux.opt.nodsatnone = coef1.aux.opt.nodsatnone;
            coef.aux.opt.gwchlint = coef1.aux.opt.gwchlint;
            coef.aux.opt.gwchnone = coef1.aux.opt.gwchnone;
            coef.aux.opt.inferaprx = coef1.aux.opt.inferaprx;
            coef.aux.opt.rmin = coef1.aux.opt.rmin;
            coef.aux.opt.method = coef1.aux.opt.method;
            coef.aux.opt.tunconst = coef1.aux.opt.tunconst;
            coef.aux.opt.tunrdn = coef1.aux.opt.tunrdn;
            coef.aux.opt.linci = coef1.aux.opt.linci;
            coef.aux.opt.white = coef1.aux.opt.white;
            coef.aux.opt.nrlzn = coef1.aux.opt.nrlzn;
            coef.aux.opt.lsfrqosmp = coef1.aux.opt.lsfrqosmp;
            coef.aux.opt.nodiagn = coef1.aux.opt.nodiagn;
            coef.aux.opt.diagnplots = coef1.aux.opt.diagnplots;
            coef.aux.opt.diagnminsnr = coef1.aux.opt.diagnminsnr;
            coef.aux.opt.ordercnstit = coef1.aux.opt.ordercnstit;
            coef.aux.opt.runtimedisp = coef1.aux.opt.runtimedisp;
            coef.aux.opt.cnstit = coef1.aux.opt.cnstit;
            if ~ndtst
                coef.diagn.name = coef1.diagn.name;
                coef.diagn.lo.name = coef1.diagn.lo.name;
                coef.diagn.lo.RR = coef1.diagn.lo.RR;
                coef.diagn.lo.RNM = coef1.diagn.lo.RNM;
                coef.diagn.hi.name = coef1.diagn.hi.name;
                coef.diagn.hi.RR = coef1.diagn.hi.RR;
                coef.diagn.hi.RNM = coef1.diagn.hi.RNM;
            end
        end
        % store results from this particular (is'th) run 
        if ~isempty(vin)            
            coef.Lsmaj(:,is) = coef1.Lsmaj;
            coef.Lsmaj_ci(:,is) = coef1.Lsmaj_ci;
            coef.Lsmin(:,is) = coef1.Lsmin;
            coef.Lsmin_ci(:,is) = coef1.Lsmin_ci;
            coef.theta(:,is) = coef1.theta;
            coef.theta_ci(:,is) = coef1.theta_ci;
            coef.umean(is) = coef1.umean;
            coef.vmean(is) = coef1.vmean;
            if ~nttst
                coef.uslope(is) = coef1.uslope;
                coef.vslope(is) = coef1.vslope;
            end 
        else
            coef.A(:,is) = coef1.A;
            coef.A_ci(:,is) = coef1.A_ci;
            coef.mean(is) = coef1.mean;
            if ~nttst
                coef.slope(is) = coef1.slope;
            end
        end
        coef.g(:,is) = coef1.g;
        coef.g_ci(:,is) = coef1.g_ci;
        coef.results{is} = coef1.results;
        coef.aux.rundescr{is} = coef1.aux.rundescr;
        coef.aux.reftime(is) = coef1.aux.reftime;
        coef.aux.opt.equi(is) = coef1.aux.opt.equi;
        if ~ndtst
            coef.diagn.PE(:,is) = coef1.diagn.PE;
            coef.diagn.SNR(:,is) = coef1.diagn.SNR;
            coef.diagn.TVraw(is) = coef1.diagn.TVraw;
            coef.diagn.TVallc(is) = coef1.diagn.TVallc;
            coef.diagn.TVsnrc(is) = coef1.diagn.TVsnrc;
            coef.diagn.PTVallc(is) = coef1.diagn.PTVallc;
            coef.diagn.PTVsnrc(is) = coef1.diagn.PTVsnrc;
            coef.diagn.SNRallc(is) = coef1.diagn.SNRallc;
            coef.diagn.K(is) = coef1.diagn.K;
            coef.diagn.lo.CorMx(:,is) = coef1.diagn.lo.CorMx;
            coef.diagn.hi.CorMx(:,is) = coef1.diagn.hi.CorMx;
            coef.diagn.table{is} = coef1.diagn.table;
        end
    end
    if ~isempty(infind)
        coef.aux.opt.infer = infer0;
    else
        coef.aux.opt.infer = [];
    end
    coef.aux.lat = lat;
    % reshape back to original shapes (from vectors) if needed
    if n_d > 1
        coef.g = reshape(coef.g,[nallc alln]);
        coef.g_ci = reshape(coef.g_ci,[nallc alln]);
        if ~isempty(vin)
            coef.Lsmaj = reshape(coef.Lsmaj,[nallc alln]);
            coef.Lsmaj_ci = reshape(coef.Lsmaj_ci,[nallc alln]);
            coef.Lsmin = reshape(coef.Lsmin,[nallc alln]);
            coef.Lsmin_ci = reshape(coef.Lsmin_ci,[nallc alln]);
            coef.theta = reshape(coef.theta,[nallc alln]);
            coef.theta_ci = reshape(coef.theta_ci,[nallc alln]);
            coef.umean = reshape(coef.umean,alln);
            coef.vmean = reshape(coef.vmean,alln);
            if ~nttst
                coef.uslope = reshape(coef.uslope,alln);
                coef.vslope = reshape(coef.vslope,alln);
            end
        else
            coef.A = reshape(coef.A,[nallc alln]);
            coef.A_ci = reshape(coef.A_ci,[nallc alln]);
            coef.mean = reshape(coef.mean,alln);
            if ~nttst
                coef.slope = reshape(coef.slope,alln);
            end
        end
        coef.results = reshape(coef.results,alln);
        if ~isempty(infind)
            if ~isempty(vin)
                if ~isequal(size(coef.aux.opt.infer.amprat),[2*ninf 1])
                    coef.aux.opt.infer.amprat = ...
                        reshape(coef.aux.opt.infer.amprat,[2*ninf alln]);
                    coef.aux.opt.infer.phsoff = ...
                        reshape(coef.aux.opt.infer.phsoff,[2*ninf alln]);
                end
            else
                if ~isequal(size(coef.aux.opt.infer.amprat),[ninf 1])
                    coef.aux.opt.infer.amprat = ...
                        reshape(coef.aux.opt.infer.amprat,[ninf alln]);
                    coef.aux.opt.infer.phsoff = ...
                        reshape(coef.aux.opt.infer.phsoff,[ninf alln]);
                end
            end
        end
        if ~isequal(size(lat),[1 1])
            coef.aux.lat = reshape(lat,alln);
        end
        coef.aux.rundescr = reshape(coef.aux.rundescr,alln);
        coef.aux.reftime = reshape(coef.aux.reftime,alln);
        coef.aux.opt.equi = reshape(coef.aux.opt.equi,alln);
        if ~ndtst
            coef.diagn.PE = reshape(coef.diagn.PE,[nallc alln]);
            coef.diagn.SNR = reshape(coef.diagn.SNR,[nallc alln]);
            coef.diagn.TVraw = reshape(coef.diagn.TVraw,alln);
            coef.diagn.TVallc = reshape(coef.diagn.TVallc,alln);
            coef.diagn.TVsnrc = reshape(coef.diagn.TVsnrc,alln);
            coef.diagn.PTVallc = reshape(coef.diagn.PTVallc,alln);
            coef.diagn.PTVsnrc = reshape(coef.diagn.PTVsnrc,alln);
            coef.diagn.SNRallc = reshape(coef.diagn.SNRallc,alln);
            coef.diagn.K = reshape(coef.diagn.K,alln);
            coef.diagn.lo.CorMx = reshape(coef.diagn.lo.CorMx,...
                [nallc alln]);
            coef.diagn.hi.CorMx = reshape(coef.diagn.hi.CorMx,...
                [nallc alln]);
            coef.diagn.table = reshape(coef.diagn.table,alln);
        end
    end
    % reorder fields
    if coef.aux.opt.twodim
        fldord = {'name'; 'Lsmaj'; 'Lsmaj_ci';...
            'Lsmin'; 'Lsmin_ci'; 'theta'; 'theta_ci'; 'g'; 'g_ci';...
            'umean'; 'vmean';};
    else
        fldord = {'name'; 'A'; 'A_ci'; 'g'; 'g_ci';...
            'mean';};
    end
    if ~coef.aux.opt.notrend
        if coef.aux.opt.twodim
            fldord{end+1} = 'uslope';
            fldord{end+1} = 'vslope';
        else
            fldord{end+1} = 'slope';
        end
    end
    fldord{end+1} = 'results';
    fldord{end+1} = 'aux';
    if ~coef.aux.opt.nodiagn
        fldord{end+1} = 'diagn';
    end
    coef = orderfields(coef,fldord);
    coef.aux = orderfields(coef.aux,{'rundescr';'opt';'frq';'lind';...
        'lat';'reftime';});    
    coef.aux.opt = orderfields(coef.aux.opt,{'twodim';'equi';'notrend';...
        'prefilt';'nodsatlint';'nodsatnone';'gwchlint';'gwchnone';...
        'infer';'inferaprx';'rmin';'method';'tunconst';'tunrdn';'linci';...
        'white';'nrlzn';'lsfrqosmp';'nodiagn';'diagnplots';...
        'diagnminsnr';'ordercnstit';'runtimedisp';'cnstit';});
    if ~ndtst
        coef.diagn = orderfields(coef.diagn,{'name';'PE';'SNR';'TVraw';...
            'TVallc';'TVsnrc';'PTVallc';'PTVsnrc';'SNRallc';'K';'lo';...
            'hi';'table';});
        coef.diagn.lo = orderfields(coef.diagn.lo,{'name';'RR';'RNM';...
            'CorMx';});
        coef.diagn.hi = orderfields(coef.diagn.hi,{'name';'RR';'RNM';...
            'CorMx';});
    end
end

%%--------------------------------------------------------- 
function coef = ut_solv1(tin,uin,vin,lat,cnstit,varargin)
% UT_SOLV1()
% Tidal analysis of a single record. See comments for UT_SOLV().
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 

fprintf('ut_solv: ');

%% input checking/conditioning, options parsing
[nt,t,u,v,tref,lor,elor,opt,tgd,uvgd] = ...
    ut_slvinit(tin,uin,vin,cnstit,varargin);
    
%% constituent selection, initial storage
opt.cnstit = cnstit;
[nNR,nR,nI,cnstit,coef] = ut_cnstitsel(tref,opt.rmin/(24*lor),...
    opt.cnstit,opt.infer);
coef.aux.rundescr = ut_rundescr(opt,nNR,nR,nI,t,tgd,uvgd,lat);
coef.aux.opt = opt;
coef.aux.lat = lat;

%% basis function matrix
fprintf('matrix prep ... ');
ngflgs = [opt.nodsatlint opt.nodsatnone opt.gwchlint opt.gwchnone];
E = ut_E(t,tref,cnstit.NR.frq,cnstit.NR.lind,lat,ngflgs,opt.prefilt);
B = [E conj(E)];
if ~isempty(opt.infer)
    Etilp = nan*ones(nt,nR);
    Etilm = Etilp;
    if ~opt.inferaprx
        for k=1:nR
            E = ut_E(t,tref,cnstit.R{k}.frq,cnstit.R{k}.lind,lat,...
                ngflgs,opt.prefilt);
            Q = ut_E(t,tref,cnstit.R{k}.I.frq,cnstit.R{k}.I.lind,lat,...
                ngflgs,opt.prefilt) ./ repmat(E,1,cnstit.R{k}.nI);
            Etilp(:,k) = E.*(1+sum(Q.*cnstit.R{k}.I.Rp(ones(nt,1),:),2));
            Etilm(:,k) = E.*(1+sum(Q.*...
                conj(cnstit.R{k}.I.Rm(ones(nt,1),:)),2));
        end
    else
        Q = nan*ones(1,nR);
        beta = Q;
        for k=1:nR
            E = ut_E(t,tref,cnstit.R{k}.frq,cnstit.R{k}.lind,lat,...
                ngflgs,opt.prefilt);
            Etilp(:,k) = E;
            Etilm(:,k) = E;
            Q(k)=ut_E(tref,tref,cnstit.R{k}.I.frq,cnstit.R{k}.I.lind,...
                lat,ngflgs,opt.prefilt)  ./ ...
                ut_E(tref,tref,cnstit.R{k}.frq,cnstit.R{k}.lind,...
                lat,ngflgs,opt.prefilt);
            arg = pi*lor*24*(cnstit.R{k}.I.frq-cnstit.R{k}.frq)*(nt+1)/nt;
            beta(k) = sin(arg)/arg;
        end    
    end
    B = [B Etilp conj(Etilm)];
end
if opt.notrend
    B = [B ones(nt,1)];
    nm = 2*(nNR + nR) + 1;
else
    B = [B ones(nt,1) (t-tref)/lor];
    nm = 2*(nNR + nR) + 2;
end

%% solution
fprintf('sol''n ... ');
xraw = u;
if opt.twodim
    xraw = complex(u,v);
end
if isequal(opt.method,'ols')
    m = B\xraw;
    W = sparse(1:nt,1:nt,1);
else
    lastwarn('');
    [m,solnstats] = robustfit(B,ctranspose(xraw),...
        opt.method,opt.tunconst,'off');
    if isequal(lastwarn,'Iteration limit reached.')
        % nan-fill, create coef.results, reorder coef fields, 
        %   do runtime display
        coef = ut_finish(coef,nNR,nR,nI,elor,cnstit);
        % abort remainder of calcs
        return;        
    end
    W = sparse(1:nt,1:nt,solnstats.w);
end
xmod = B*m;
if ~opt.twodim
    xmod = real(xmod);
end
e = W*(xraw-xmod);

%% NR & R cnstits: complex & cos/sin coeffs; current ellipse params
nc = nNR+nR;
ap = m([1:nNR 2*nNR+(1:nR)]);
am = m([nNR+(1:nNR) 2*nNR+nR+(1:nR)]);
if ~isempty(opt.infer)
    if opt.inferaprx
        for k=1:nR
              ap(nNR+k) = ap(nNR+k)./(1+beta(k)*cnstit.R{k}.I.Rp.*Q(k));
              am(nNR+k) = am(nNR+k)./...
                  (1+beta(k)*cnstit.R{k}.I.Rm.*conj(Q(k)));
        end
    end
end
Xu = real(ap+am);
Yu = -imag(ap-am);
if ~opt.twodim
    [coef.A,~,~,coef.g] = ut_cs2cep([Xu Yu]);
else
    Xv = imag(ap+am);
    Yv = real(ap-am);
    [coef.Lsmaj,coef.Lsmin,coef.theta,coef.g] = ut_cs2cep([Xu Yu Xv Yv]);
end
    
%% mean and trend
if opt.twodim
    if opt.notrend
        coef.umean = real(m(end));
        coef.vmean = imag(m(end));
    else
        coef.umean = real(m(end-1));
        coef.vmean = imag(m(end-1));
        coef.uslope = real(m(end))/lor;
        coef.vslope = imag(m(end))/lor;
    end
else
    if opt.notrend
        coef.mean = real(m(end));
    else
        coef.mean = real(m(end-1));
        coef.slope = real(m(end))/lor;
    end
end

%% spectral power (Puu, Pvv, Puv) of residual
fprintf('conf. int''vls ... ');
if ~opt.white
    % band-averaged (ba) spectral densities
    if opt.equi
        if sum(tgd)>sum(uvgd)  
            efill = interp1(t,e,tin(tgd));
            if any(isnan(efill)) % fill start&/end nans w/ nearest good
                ind = find(isnan(efill));
                ind2 = ind(ind<find(~isnan(efill),1,'first')); 
                efill(ind2) = efill(max(ind2)+1);
                ind2 = ind(ind>find(~isnan(efill),1,'last'));
                efill(ind2) = efill(min(ind2)-1);
            end
            ba = ut_pdgm(tin(tgd),efill,coef.aux.frq,1,0);
        else
            ba = ut_pdgm(tin(tgd),e,coef.aux.frq,1,0);
        end
    else
        ba = ut_pdgm(t,e,coef.aux.frq,0,opt.lsfrqosmp);
    end
    % power [ (e units)^2 ] from spectral density [ (e units)^2 / cph ]
    df = 1/(elor*24);
    ba.Puu = ba.Puu*df;
    if opt.twodim
        ba.Pvv = ba.Pvv*df;
        ba.Puv = ba.Puv*df;
    end
    % assign band-avg power values to NR & R freqs
    Puu = zeros(size(coef.aux.frq));
    if opt.twodim
        Pvv = Puu;
        Puv = Pvv;
    end
    for i = 1:length(ba.Puu)
        ind = find(coef.aux.frq>=ba.fbnd(i,1) & ...
            coef.aux.frq<=ba.fbnd(i,2));
        Puu(ind) = ba.Puu(i);
        if opt.twodim
            Pvv(ind) = ba.Pvv(i);
            Puv(ind) = ba.Puv(i);
        end
    end
end

%% NR & R cnstits: confidence intervals
varMSM = real((ctranspose(xraw)*W*xraw - ...
    ctranspose(m)*ctranspose(B)*W*xraw)/(nt-nm));
gamC = inv(ctranspose(B)*W*B)*varMSM; %#ok
gamP = inv(transpose(B)*W*B)*((transpose(xraw)*W*xraw - ...
    transpose(m)*transpose(B)*W*xraw)/(nt-nm)); %#ok
Gall = gamC+gamP;
Hall = gamC-gamP;
% inits
coef.g_ci = nan*ones(size(coef.g));
if opt.twodim
    coef.Lsmaj_ci = coef.g_ci;
    coef.Lsmin_ci = coef.g_ci;
    coef.theta_ci = coef.g_ci;
    varcov_mCw = nan*ones(nc,4,4);
else
    coef.A_ci = coef.g_ci;
    varcov_mCw = nan*ones(nc,2,2);
end
if ~opt.white
    varcov_mCc = varcov_mCw;
end
% main loop
for c=1:nc 
    G = [Gall(c,c) Gall(c,c+nc); Gall(c+nc,c) Gall(c+nc,c+nc);];
    H = [Hall(c,c) Hall(c,c+nc); Hall(c+nc,c) Hall(c+nc,c+nc);];
    varXu = real(G(1,1)+G(2,2)+2*G(1,2))/2;
    varYu = real(H(1,1)+H(2,2)-2*H(1,2))/2;
    if opt.twodim
        varXv = real(H(1,1)+H(2,2)+2*H(1,2))/2;
        varYv = real(G(1,1)+G(2,2)-2*G(1,2))/2;
    end
    if opt.linci % linearized
        if ~opt.twodim
            varcov_mCw(c,:,:) = diag([varXu varYu]);
            if ~opt.white
                den = varXu + varYu;
                varXu = Puu(c)*varXu/den;
                varYu = Puu(c)*varYu/den;
                varcov_mCc(c,:,:) = diag([varXu varYu]);
            end
            [sig1,sig2]= ut_linci(Xu(c),Yu(c),sqrt(varXu),sqrt(varYu));
            coef.A_ci(c) = 1.96*sig1;
            coef.g_ci(c) = 1.96*sig2;
        else
            varcov_mCw(c,:,:) = diag([varXu varYu varXv varYv]);
            if ~opt.white
                den = varXv + varYv;
                varXv = Pvv(c)*varXv/den;
                varYv = Pvv(c)*varYv/den;
                varcov_mCc(c,:,:) = diag([varXu varYu varXv varYv]);
            end
            [sig1,sig2] = ut_linci(Xu(c)+1i*Xv(c),Yu(c)+1i*Yv(c),...
                sqrt(varXu)+1i*sqrt(varXv),sqrt(varYu)+1i*sqrt(varYv));
            coef.Lsmaj_ci(c) = 1.96*real(sig1);
            coef.Lsmin_ci(c) = 1.96*imag(sig1);
            coef.g_ci(c) = 1.96*real(sig2);
            coef.theta_ci(c) = 1.96*imag(sig2);
        end
    else % monte carlo
        covXuYu = imag(H(1,1)-H(1,2)+H(2,1)-H(2,2))/2;
        Duu = [varXu covXuYu; covXuYu varYu;];
        varcov_mCw(c,1:2,1:2) = Duu;
        if ~opt.white
            Duu = Puu(c)*Duu/trace(Duu);
            varcov_mCc(c,1:2,1:2) = Duu;
        end
        if ~opt.twodim
            if ~opt.white
                varcov_mCc(c,:,:) = ut_nearposdef(squeeze(...
                    varcov_mCc(c,:,:)));
                mCall = mvnrnd([Xu(c) Yu(c)],...
                    squeeze(varcov_mCc(c,:,:)),opt.nrlzn);
            else
                mCall = mvnrnd([Xu(c) Yu(c)],squeeze(varcov_mCw(c,:,:)),...
                    opt.nrlzn);    
            end
            [A,~,~,g] = ut_cs2cep(mCall);
            coef.A_ci(c) = 1.96*median(abs(A-median(A)))/0.6745;
            g(1) = coef.g(c);
            g = ut_cluster(g,360);
            coef.g_ci(c) = 1.96*median(abs(g-median(g)))/0.6745;
        else
            covXvYv = imag(G(1,1)-G(1,2)+G(2,1)-G(2,2))/2;
            Dvv = [varXv covXvYv; covXvYv varYv;];
            varcov_mCw(c,3:4,3:4) = Dvv;
            if ~opt.white
                Dvv = Pvv(c)*Dvv/trace(Dvv);
                varcov_mCc(c,3:4,3:4) = Dvv;
            end
            covXuXv = imag(-H(1,1)-H(1,2)-H(2,1)-H(2,2))/2;
            covXuYv = real(G(1,1)-G(2,2))/2;
            covYuXv = real(-H(1,1)+H(2,2))/2;
            covYuYv = imag(-G(1,1)+G(1,2)+G(2,1)-G(2,2))/2;
            Duv = [covXuXv covXuYv; covYuXv covYuYv;];
            varcov_mCw(c,1:2,3:4) = Duv;
            varcov_mCw(c,3:4,1:2) = transpose(Duv);
            if ~opt.white
                if sum(abs(Duv(:)))>0
                    Duv = Puv(c)*Duv/sum(abs(Duv(:)));
                    varcov_mCc(c,1:2,3:4) = Duv;
                    varcov_mCc(c,3:4,1:2) = transpose(Duv);
                else
                    varcov_mCc(c,1:2,3:4) = 0;
                    varcov_mCc(c,3:4,1:2) = 0;
                end
                varcov_mCc(c,:,:) = ut_nearposdef(squeeze(...
                    varcov_mCc(c,:,:)));
                mCall = mvnrnd([Xu(c) Yu(c) Xv(c) Yv(c)],...
                    squeeze(varcov_mCc(c,:,:)),opt.nrlzn);
            else
                mCall = mvnrnd([Xu(c) Yu(c) Xv(c) Yv(c)],...
                    squeeze(varcov_mCw(c,:,:)),opt.nrlzn);
            end
            [Lsmaj,Lsmin,theta,g]= ut_cs2cep(mCall);
            coef.Lsmaj_ci(c) =1.96*median(abs(Lsmaj-median(Lsmaj)))/0.6745;
            coef.Lsmin_ci(c) =1.96*median(abs(Lsmin-median(Lsmin)))/0.6745;
            theta(1) = coef.theta(c);
            theta = ut_cluster(theta,360);
            coef.theta_ci(c) =1.96*median(abs(theta-median(theta)))/0.6745;
            g(1) = coef.g(c);
            g = ut_cluster(g,360);
            coef.g_ci(c) = 1.96*median(abs(g-median(g)))/0.6745;
        end
    end
end

%% I cnstits: complex & cos/sin coeffs, c.e.p., & conf ints
if ~isempty(opt.infer)
    % complex coeffs
    apI = nan*ones(1,nI);
    amI = apI;
    ind = 1;
    for k=1:nR
        apI(ind:ind+cnstit.R{k}.nI-1) = cnstit.R{k}.I.Rp*ap(nNR+k);
        amI(ind:ind+cnstit.R{k}.nI-1) = cnstit.R{k}.I.Rm*am(nNR+k);
        ind = ind + cnstit.R{k}.nI;
    end    
    % cos/sin coeffs and c.e.p.
    XuI = real(apI+amI)';
    YuI = -imag(apI-amI)';
    if ~opt.twodim
        [A,~,~,g] = ut_cs2cep([XuI YuI]);
        coef.A = [coef.A; A;];
        coef.g = [coef.g; g;];
    else
        XvI = imag(apI+amI)';
        YvI = real(apI-amI)';
        [Lsmaj,Lsmin,theta,g] = ut_cs2cep([XuI YuI XvI YvI]);
        coef.Lsmaj = [coef.Lsmaj; Lsmaj;];
        coef.Lsmin = [coef.Lsmin; Lsmin;];
        coef.theta = [coef.theta; theta;];
        coef.g = [coef.g; g;];
    end
    % conf ints
    if opt.linci % linearized
        for k=1:nR
            if opt.white
                if ~opt.twodim
                    varReap = 0.25*varcov_mCw(nNR+k,1,1);
                    varImap = 0.25*varcov_mCw(nNR+k,2,2);
                else
                    varReap = 0.25*(varcov_mCw(nNR+k,1,1) + ...
                        varcov_mCw(nNR+k,4,4));
                    varImap = 0.25*(varcov_mCw(nNR+k,2,2) + ...
                        varcov_mCw(nNR+k,3,3));
                end
            else
                if ~opt.twodim
                    varReap = 0.25*varcov_mCc(nNR+k,1,1);
                    varImap = 0.25*varcov_mCc(nNR+k,2,2);
                else
                    varReap = 0.25*(varcov_mCc(nNR+k,1,1) + ...
                        varcov_mCc(nNR+k,4,4));
                    varImap = 0.25*(varcov_mCc(nNR+k,2,2) + ...
                        varcov_mCc(nNR+k,3,3));
                end
            end
            varXuHH = (real(cnstit.R{k}.I.Rp).^2 + ...
                real(cnstit.R{k}.I.Rm).^2)*varReap + ...
                (imag(cnstit.R{k}.I.Rp).^2 + ...
                imag(cnstit.R{k}.I.Rm).^2)*varImap;
            varYuHH = (imag(cnstit.R{k}.I.Rp).^2 + ...
                imag(cnstit.R{k}.I.Rm).^2)*varReap + ...
                (real(cnstit.R{k}.I.Rp).^2 + ...
                real(cnstit.R{k}.I.Rm).^2)*varImap;
            for c = 1:length(varXuHH)
                if ~opt.twodim
                    [sig1,sig2]= ut_linci(Xu(nNR+k),Yu(nNR+k),...
                        sqrt(varXuHH(c)),sqrt(varYuHH(c)));
                    coef.A_ci = [coef.A_ci; 1.96*sig1;];
                    coef.g_ci = [coef.g_ci; 1.96*sig2;];
                else
                    [sig1,sig2]= ut_linci(Xu(nNR+k)+1i*Xv(nNR+k),...
                        Yu(nNR+k)+1i*Yv(nNR+k),...
                        sqrt(varXuHH(c))+1i*sqrt(varYuHH(c)),...
                        sqrt(varYuHH(c))+1i*sqrt(varXuHH(c)));
                    coef.Lsmaj_ci = [coef.Lsmaj_ci; 1.96*real(sig1);];
                    coef.Lsmin_ci = [coef.Lsmin_ci; 1.96*imag(sig1);];
                    coef.g_ci = [coef.g_ci; 1.96*real(sig2);];
                    coef.theta_ci = [coef.theta_ci; 1.96*imag(sig2);];
                end
            end
        end
    else % monte carlo
        ind = 0;
        for k = 1:nR
            if ~opt.twodim
                if opt.white
                    mCall = mvnrnd([Xu(nNR+k) Yu(nNR+k)],...
                        squeeze(varcov_mCw(nNR+k,:,:)),opt.nrlzn);
                else
                    mCall = mvnrnd([Xu(nNR+k) Yu(nNR+k)],...
                        squeeze(varcov_mCc(nNR+k,:,:)),opt.nrlzn);
                end
                [A,~,~,g] = ut_cs2cep(mCall);
                g(1) = coef.g(nNR+k);
                g = ut_cluster(g,360);
                for lk = 1:cnstit.R{k}.nI
                    ind = ind+1;
                    apI = 0.5*cnstit.R{k}.I.Rp(lk)*(A.*exp(-1i*g*pi/180));
                    amI = conj(apI);
                    XuI = real(apI+amI);
                    YuI = -imag(apI-amI);
                    [A,~,~,g] = ut_cs2cep([XuI YuI]);
                    A_ci = median(abs(A-median(A)))/0.6745;
                    coef.A_ci = [coef.A_ci; 1.96*A_ci;];
                    g = ut_cluster(g,360);
                    g_ci = median(abs(g-median(g)))/0.6745;
                    coef.g_ci = [coef.g_ci; 1.96*g_ci;];
                end
            else
                if opt.white
                    mCall = mvnrnd([Xu(nNR+k) Yu(nNR+k) Xv(nNR+k) ...
                        Yv(nNR+k)],squeeze(varcov_mCw(nNR+k,:,:)),...
                        opt.nrlzn);
                else
                    mCall = mvnrnd([Xu(nNR+k) Yu(nNR+k) Xv(nNR+k) ...
                        Yv(nNR+k)],squeeze(varcov_mCc(nNR+k,:,:)),...
                        opt.nrlzn);
                end                    
                [Lsmaj,Lsmin,theta,g] = ut_cs2cep(mCall);
                g(1) = coef.g(nNR+k);
                g = ut_cluster(g,360);
                theta(1) = coef.theta(nNR+k);
                theta = ut_cluster(theta,360);
                for lk = 1:cnstit.R{k}.nI
                    ind = ind+1;
                    apI = cnstit.R{k}.I.Rp(lk).*0.5*(Lsmaj+Lsmin).* ...
                        exp(1i*(theta-g)*pi/180);
                    amI = cnstit.R{k}.I.Rm(lk).*0.5*(Lsmaj-Lsmin).* ...
                        exp(1i*(theta+g)*pi/180);
                    XuI = real(apI+amI);
                    YuI = -imag(apI-amI);
                    XvI = imag(apI+amI);
                    YvI = real(apI-amI);
                    [Lsmaj,Lsmin,theta,g] = ut_cs2cep([XuI YuI XvI YvI]);
                    Lsmaj_ci = median(abs(Lsmaj-median(Lsmaj)))/0.6745;
                    coef.Lsmaj_ci = [coef.Lsmaj_ci; 1.96*Lsmaj_ci;];
                    Lsmin_ci = median(abs(Lsmin-median(Lsmin)))/0.6745;
                    coef.Lsmin_ci = [coef.Lsmin_ci; 1.96*Lsmin_ci;];
                    theta = ut_cluster(theta,360);
                    theta_ci = median(abs(theta-median(theta)))/0.6745;
                    coef.theta_ci = [coef.theta_ci; 1.96*theta_ci;];
                    g = ut_cluster(g,360);
                    g_ci = median(abs(g-median(g)))/0.6745;
                    coef.g_ci = [coef.g_ci; 1.96*g_ci;];
                end
            end            
        end
    end
end

%% diagnostics
if ~opt.nodiagn
    fprintf('diagnostics ... ');
    if opt.twodim
        PE = sum(coef.Lsmaj.^2 + coef.Lsmin.^2);
        PE = 100*(coef.Lsmaj.^2 + coef.Lsmin.^2)/PE;
        SNR = (coef.Lsmaj.^2 +coef.Lsmin.^2)./...
            ((coef.Lsmaj_ci/1.96).^2 + (coef.Lsmin_ci/1.96).^2);
    else
        PE = 100*coef.A.^2/sum(coef.A.^2);
        SNR = (coef.A.^2)./((coef.A_ci/1.96).^2);
    end
    [~,indPE] = sort(PE,'descend');
    coef.diagn.name = coef.name(indPE);
    coef.diagn.PE = PE(indPE);
    coef.diagn.SNR = SNR; % used in ut_diagntable; ordered by PE there
    if opt.twodim
        [coef.diagn,usnrc,vsnrc] = ut_diagntable(coef,cnstit,...
            t,u,v,xmod,m,B,W,varMSM,Gall,Hall,elor,varcov_mCw,indPE);
    else
        [coef.diagn,usnrc,~] = ut_diagntable(coef,cnstit,...
            t,u,[],xmod,m,B,W,varMSM,Gall,Hall,elor,varcov_mCw,indPE);
    end
    if opt.diagnplots
        tmp = nan*ones(size(uin));
        tmp(uvgd) = usnrc;
        usnrc = tmp;
        tmp = nan*ones(size(uin));
        tmp(uvgd) = e;
        e = tmp;
        if opt.twodim
            tmp = nan*ones(size(uin));
            tmp(uvgd) = vsnrc;
            vsnrc = tmp;
            ut_diagnfigs(coef,indPE,tin,uin,vin,usnrc,vsnrc,e);
        else
            ut_diagnfigs(coef,indPE,tin,uin,[],usnrc,[],e);
        end
    end
end

%% re-order constituents
if ~isempty(opt.ordercnstit)
    if isequal(opt.ordercnstit,'frq')
        [~,ind] = sort(coef.aux.frq);
    elseif isequal(opt.ordercnstit,'snr')
        if ~opt.nodiagn
            [~,ind] = sort(coef.diagn.SNR,'descend');
        else
            if opt.twodim
                SNR = (coef.Lsmaj.^2 +coef.Lsmin.^2)./...
                    ((coef.Lsmaj_ci/1.96).^2 + (coef.Lsmin_ci/1.96).^2);
            else
                SNR = (coef.A.^2)./((coef.A_ci/1.96).^2);
            end
            [~,ind] = sort(SNR,'descend');
        end
    else        
        [~,ind] = ismember(cellstr(opt.ordercnstit),cellstr(opt.cnstit));
    end
else
    if ~opt.nodiagn
        ind = indPE;
    else
        if opt.twodim
            PE = sum(coef.Lsmaj.^2 + coef.Lsmin.^2);
            PE = 100*(coef.Lsmaj.^2 + coef.Lsmin.^2)/PE;
        else
            PE = 100*coef.A.^2/sum(coef.A.^2);
        end
        [~,ind] = sort(PE,'descend');
    end
end
coef.name = coef.name(ind);
coef.g = coef.g(ind);
coef.g_ci = coef.g_ci(ind);
if opt.twodim
    coef.Lsmaj = coef.Lsmaj(ind);
    coef.Lsmaj_ci = coef.Lsmaj_ci(ind);
    coef.Lsmin = coef.Lsmin(ind);
    coef.Lsmin_ci = coef.Lsmin_ci(ind);
    coef.theta = coef.theta(ind);
    coef.theta_ci = coef.theta_ci(ind);
else
    coef.A = coef.A(ind);
    coef.A_ci = coef.A_ci(ind);
end
coef.aux.frq = coef.aux.frq(ind);
coef.aux.lind = coef.aux.lind(ind);
fprintf('done.\n');

%% create coef.results, reorder coef fields, do runtime display
coef = ut_finish(coef,nNR,nR,nI,elor,cnstit);

%%--------------------------------------------------------- 
function [nt,t,u,v,tref,lor,elor,opt,tgd,uvgd] = ...
    ut_slvinit(tin,uin,vin,cnstit,args)
% UT_SLVINIT()
% initial input parsing and conditioning for UT_SOLV1()
% inputs
%   tin = input times [datenum UTC] passed in to ut_solv (ntin x 1)
%   uin, vin = raw inputs [units arb.] passed in to ut_solv (ntin x 1)
%   cnstit = cell array of 4-char strings passed in to ut_solv (nc x 1)
%   args = varargin from call to ut_solv
% outputs
%   nt = number of t/u/v values ( b/c of nans, can be < length(tin) )
%   t = tin values having non-nan uin, vin [datenum UTC] (nt x 1)
%   u,v = uin,vin elements [units arb.] where both are non-nan (nt x 1)
%   tref = central reference time [datenum UTC] (1 x 1)
%   lor = length of record, (max(t) - min(t)) [days] (1 x 1)
%   elor = effective length of record, lor*nt/(nt-1) [days] (1 x 1)
%   opt = option flags information
%   tgd = result of ~isnan(tin) (ntin x 1)
%   uvgd = result of ~isnan(uin) & tgd [ & ~isnan(vin) ] (ntin x 1)
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 

%% input checking/reshaping
if ~isequal(size(tin),size(uin))
    error(['ut_solv: vectors of input times and input values '...
        'must be same size.']);
end
tgd = ~isnan(tin);
uin = uin(tgd);
if ~isempty(vin)
    if ~isequal(size(tin),size(vin))
        error(['ut_solv: vectors of input times and input values '...
            'must be same size.']);
    end
    vin = vin(tgd);
end
tin = tin(tgd);
tin = tin(:);
if ~isequal(sort(tin),tin) || ~isequal(length(unique(tin)),length(tin))
    error('ut_solv: input times must be monotonically increasing.');
end
if min(tin)<datenum(500,0,0) || max(tin)>datenum(3500,0,0)
    error('ut_solv: input times should be from the present millenium.');
end
uin = uin(:);
if ~isequal(size(uin,2),1) || size(uin,1)<1
    error('ut_solv: multiple input values required.');
end
if isempty(vin)
    opt.twodim = 0;
    v = [];
else
    opt.twodim = 1;
    vin = vin(:);
end
% only keep t,u,v values with non-nan u(&v) 
if opt.twodim
    uvgd = ~isnan(uin) & ~isnan(vin);
    v = vin(uvgd);
else
    uvgd = ~isnan(uin);
end
t = tin(uvgd);
nt = length(t);
u = uin(uvgd);
if var(unique(diff(tin)))<eps
    opt.equi = 1; % based on times; u/v can still have nans ("gappy")
    lor = (max(tin)-min(tin));
    elor = lor*length(tin)/(length(tin)-1);
    tref = 0.5*(tin(1)+tin(end));
else
    opt.equi = 0;
    lor = (max(t) - min(t));
    elor = lor*nt/(nt-1);
    tref = 0.5*(t(1)+t(end));
end
%% options
opt.notrend = 0;
opt.prefilt = [];
opt.nodsatlint = 0;
opt.nodsatnone = 0;
opt.gwchlint = 0;
opt.gwchnone = 0;
opt.infer = [];
opt.inferaprx = 0;
opt.rmin = 1;
opt.method = 'cauchy';
opt.tunrdn = 1;
opt.linci = 0;
opt.white = 0;
opt.nrlzn = 200; 
opt.lsfrqosmp = 1;
opt.nodiagn = 0;
opt.diagnplots = 0;
opt.diagnminsnr = 2;
opt.ordercnstit = [];
opt.runtimedisp = 'yyy';
methnotset = 1;
allmethods = {'ols','andrews','bisquare','fair','huber',...
                'logistic','talwar','welsch'};
while ~isempty(args)
    switch(lower(args{1}))
        case 'notrend'
            opt.notrend = 1;
            args(1) = [];
        case 'prefilt'
            opt.prefilt = args{2};
            args(1:2) = [];
        case 'nodsatlint'
            opt.nodsatlint = 1;
            args(1) = [];
        case 'nodsatnone'
            opt.nodsatnone = 1;
            args(1) = [];
        case 'gwchlint'
            opt.gwchlint = 1;
            args(1) = [];
        case 'gwchnone'
            opt.gwchnone = 1;
            args(1) = [];
        case 'infer'
            opt.infer = args{2};
            args(1:2) = [];
        case 'inferaprx'
            opt.inferaprx = 1;
            args(1) = [];    
        case 'rmin'
            opt.rmin = args{2};
            args(1:2) = [];
        case allmethods
            if methnotset
                opt.method = lower(args{1});
                args(1) = [];
                methnotset = 0;
            else
                error('ut_solv: only one ''method'' option allowed.');
            end
        case 'tunrdn'
            opt.tunrdn = args{2};
            args(1:2) = [];
        case 'linci'
            opt.linci = 1;
            args(1) = [];
        case 'white'
            opt.white = 1;
            args(1) = [];
        case 'nrlzn'
            opt.nrlzn = args{2};
            args(1:2) = [];
        case 'lsfrqosmp'
            opt.lsfrqosmp = round(args{2});
            args(1:2) = [];
        case 'nodiagn'
            opt.nodiagn = 1;
            args(1) = [];
        case 'diagnplots'
            opt.diagnplots = 1;
            args(1) = [];
        case 'diagnminsnr'
            opt.diagnminsnr = args{2};
            args(1:2) = [];
        case 'ordercnstit'
            opt.ordercnstit = args{2};
            args(1:2) = [];
        case 'runtimedisp'
            opt.runtimedisp = lower(args{2});
            args(1:2) = [];
        otherwise
            error(['ut_solv: unrecognized input: ' args{1}]);
    end
end
if (opt.nodsatlint + opt.nodsatnone)>1
    error(['ut_solv: only one of ''NodSatNone'' or ''NodSatLinT'' is'...
        ' allowed.']);
end
if (opt.gwchlint + opt.gwchnone)>1
    error('ut_solv: inconsistent inputs (only one gwch option allowed).');
end
if opt.lsfrqosmp<1
    error(['ut_solv: Lomb-Scargle frequency oversampling factor must be'...
        ' >=1.']);
end
if ~isempty(opt.ordercnstit) 
    if ~isequal(opt.ordercnstit,'frq') && ~isequal(opt.ordercnstit,'snr') 
        if isequal(lower(cnstit),'auto')
             error(['ut_solv: OrderCnstit argument must be ''frq'' '...
                'or ''snr'' since cnstit=''auto''.']);
        end
        if ~isequal(sort(opt.ordercnstit),sort(cnstit))
            error(['ut_solv: OrderCnstit argument (if not ''SNR'' '...
                'nor ''Frq'') must include same elements as ''cnstit''.']);
        end
    end
end
if opt.nodiagn && isequal(opt.ordercnstit,'snr')
    error(['ut_solv: OrderCnstit cannot be ''snr'' if NoDiagn is '...
        'selected.']);
end
if ~isempty(opt.infer)
    ninf = length(opt.infer.infnam);
    if opt.twodim 
        if ~isequal(size(opt.infer.amprat,1),2*ninf)...
                || ~isequal(size(opt.infer.amprat,1),2*ninf)
            error(['ut_solv: opt.infer.amprat/phsoff must each'...
                ' have 2*length(opt.infer.infnam) elements.']);
        end
    else
        if ~isequal(size(opt.infer.amprat,1),ninf)...
                || ~isequal(size(opt.infer.amprat,1),ninf)
            error(['ut_solv: opt.infer.amprat/phsoff must each'...
                ' have length(opt.infer.infnam) elements.']);
        end
    end
end
if ~isempty(opt.infer) && opt.inferaprx
    if size(opt.infer.refnam,1)>length(unique(cellstr(opt.infer.refnam)))
        error(['ut_solv: cannot infer multiple constituents from one '...
            'reference constituent with ''InfrAprx''.']);
    end
    if ~(opt.nodsatlint||opt.nodsatnone)||~(opt.gwchlint||opt.gwchnone)
        error(['ut_solv: ''InfrAprx'' requires ''nodsatlint'' or '...
            '''nodsatnone'', as well as ''gwchlint'' or ''gwchnone''.']);
    end
end
if ~isequal(opt.method,'cauchy')
    [~,ind] = ismember(opt.method,allmethods);
    allconst = [NaN 1.339 4.685 1.400 1.345 1.205 2.795 2.985];
    opt.tunconst = allconst(ind);
else
    opt.tunconst = 2.385;
end
opt.tunconst = opt.tunconst/opt.tunrdn;
if (opt.nodiagn + opt.diagnplots)>1
    error(['ut_solv: cannot use both ''NoDiagn'' and ''DiagnPlots'' '...
        'options.']);
end
nf = length(fieldnames(opt));
opt = orderfields(opt,[1:12 nf 13:(nf-1)]);
if ~isequal(length(opt.runtimedisp),3)
    error('ut_solv: ''RunTimeDisp'' must be a 3-character string.');
else
    if ~ismember(opt.runtimedisp(1),{'y','n'}) || ...
            ~ismember(opt.runtimedisp(2),{'y','n'}) || ...
            ~ismember(opt.runtimedisp(3),{'y','n'})
        error(['ut_solv: ''RunTimeDisp'', while not case-sensitive, '...
            'must include only y or n characters.']);
    end
end

%%--------------------------------------------------------- 
function [nNR,nR,nI,cnstit,coef] = ut_cnstitsel(tref,minres,incnstit,infer)
% UT_CNSTITSEL()
% carry out constituent selection
% inputs
%   tref = reference time (datenum UTC)
%   minres = freq separation (cph) used in decision tree
%   incnstit = 'cnstit' input to ut_solv
%   infer = 'opt.infer' input to ut_solv
% outputs
%   nNR,nR,nI = number non-reference, reference, inferred constituents
%   cnstit.NR.name = cellstr of 4-char names of NR constits
%   cnstit.NR.frq = frequencies (cph) of NR constits
%   cnstit.NR.lind = list indices (in ut_constants.mat) of NR constits
%   cnstit.R = empty if no inference; otherwise, for each (i'th) R constit:
%       cnstit.R{i}.name, .frq, .lind = as above, but for R constits
%       cnstit.R{i}.I{j}.name, .frq, .lind = as above for j'th I constit
%   coef.name = cellstr of names of all constituents (NR, R, and I)
%   coef.aux.frq = frequencies (cph) of all constituents
%   coef.aux.lind = list indices of all constituents
%   coef.aux.reftime = tref
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 

%% freq vals, all avail cnstits
load('ut_constants.mat','const','shallow');
[~,ader] = ut_astron(tref);
ii = isfinite(const.ishallow); %#ok
const.freq(~ii) = const.doodson(~ii,:)*ader/24;
for k = find(ii)'
    ik = const.ishallow(k) + (0:const.nshallow(k)-1);
    const.freq(k)=sum(const.freq(shallow.iname(ik)).*shallow.coef(ik));
end
%% cnstit.NR
if isequal(lower(incnstit),'auto')
    cnstit.NR.lind = find(const.df>=minres);
else
    cnstit.NR.lind = nan*ones(size(incnstit,1),1);
    for j=1:length(incnstit)
        lind1=strmatch(incnstit{j},const.name);
        if isempty(lind1),
            error(['ut_solv: unrecognized non-reference constituent: '...
                incnstit(j,:) '.']);
        elseif lind1==1,
            error(['ut_solv: do not pass in constituent Z0, '...
                'mean is removed by default.']);
        else
            cnstit.NR.lind(j) = lind1;
        end
    end
    [~,seq]=sort(const.freq(cnstit.NR.lind));
    cnstit.NR.lind = cnstit.NR.lind(seq);
end
if isempty(cnstit.NR.lind)
    error('ut_solv: no constituents specified nor auto-selected.');
end
if ~isempty(infer)
    % remove from NR any R or I constits
    nam = cellstr(char(unique([infer.infnam; infer.refnam;])));
    [tst,ind] = ismember(nam,cellstr(const.name(cnstit.NR.lind,:)));
    tst = find(tst);
    if ~isempty(tst)
        cnstit.NR.lind(ind(tst)) = [];
    end
end
cnstit.NR.frq = const.freq(cnstit.NR.lind);
cnstit.NR.name = cellstr(const.name(cnstit.NR.lind,:));
nNR = length(cnstit.NR.frq); 
%% cnstit.R 
nR = 0;
nI = 0;
cnstit.R = [];
if ~isempty(infer)
    nI = length(infer.infnam);
    allrefnames = unique(infer.refnam);
    nR = length(allrefnames);
    for k=1:nR
        % reference
        cnstit.R{k}.name = allrefnames{k};
        lind1=strmatch(cnstit.R{k}.name,const.name);
        if isempty(lind1),
            error(['ut_solv: unrecognized reference constituent: '...
                cnstit.R{k}.name '.']);
        elseif lind1==1,
            error(['ut_solv: do not pass in constituent Z0, '...
                'mean removal is handled automatically.']);
        else
            cnstit.R{k}.lind = lind1;
        end
        cnstit.R{k}.frq = const.freq(cnstit.R{k}.lind);
        % inferred
        ind = strmatch(cnstit.R{k}.name,infer.refnam);
        cnstit.R{k}.nI = length(ind);
        for lk = 1:cnstit.R{k}.nI
            cnstit.R{k}.I.Rp(lk) = infer.amprat(ind(lk)) .* ...
                exp( 1i * infer.phsoff(ind(lk))*pi/180 );
            if length(infer.amprat)>size(infer.infnam,1)
                cnstit.R{k}.I.Rm(lk)= infer.amprat(ind(lk)+nI).*...
                    exp( -1i * infer.phsoff(ind(lk)+nI)*pi/180 );
            else
                cnstit.R{k}.I.Rm(lk) = conj(cnstit.R{k}.I.Rp(lk));
            end
            cnstit.R{k}.I.name{lk} = infer.infnam{ind(lk)};
            lind1=strmatch(cnstit.R{k}.I.name{lk},const.name);
            if isempty(lind1)
                error(['ut_solv: unrecognized inference constituent: '...
                    cnstit.R{k}.I.name{lk} '.']);
            elseif lind1==1
                error(['ut_solv: do not pass in constituent Z0, '...
                    'mean removal is handled automatically.']);
            else
                cnstit.R{k}.I.lind(lk) = lind1;
            end
            cnstit.R{k}.I.frq(lk) = const.freq(cnstit.R{k}.I.lind(lk));
        end
    end
end
nallc = nNR+nR+nI;
coef.name = cnstit.NR.name; 
coef.aux.frq = nan*ones(nallc,1);
coef.aux.lind = coef.aux.frq;
coef.aux.frq(1:nNR) = cnstit.NR.frq;
coef.aux.lind(1:nNR) = cnstit.NR.lind;
if ~isempty(infer)
    for k=1:nR
        coef.aux.frq(nNR+k) = cnstit.R{k}.frq; 
        coef.aux.lind(nNR+k) = cnstit.R{k}.lind;
        coef.name(nNR+k) = cellstr(cnstit.R{k}.name); 
    end
    ind = nNR+nR+1;
    for k=1:nR
        coef.aux.frq(ind:ind+cnstit.R{k}.nI-1) = cnstit.R{k}.I.frq; 
        coef.aux.lind(ind:ind+cnstit.R{k}.nI-1) = cnstit.R{k}.I.lind; 
        coef.name(ind:ind+cnstit.R{k}.nI-1) = ...
            cellstr(char(cnstit.R{k}.I.name));
        ind = ind + cnstit.R{k}.nI;
    end
end
coef.aux.reftime = tref;

%%--------------------------------------------------------- 
function rundescr = ut_rundescr(opt,nNR,nR,nI,t,tgd,uvgd,lat)
% UT_RUNDESCR()
% build character array with descriptions of run
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 

rd = cell(10+nI,1);
rd{1} = sprintf('UTide Run Description (executed %s):',datestr(now));
if opt.equi
    if sum(tgd)>sum(uvgd)
        rd{2} = sprintf(['  Input(%dD): %d values (gappy); %d '...
            'equisp. times, dt=%#.2g hr.'],...
            opt.twodim+1,sum(uvgd),sum(tgd),24*mode(diff(t)));
    else
        rd{2} = sprintf(['  Input(%dD): %d values, equispaced times,'...
            ' dt=%#.3g hr.'],...
            opt.twodim+1,sum(uvgd),24*mode(diff(t)));
    end
else
    rd{2} = sprintf(['  Input(%dD): %d vals; irreg., dt(min/mean/max)='...
        ' %#.2g/%#.2g/%#.2g hr.'],opt.twodim+1,length(t),...
        24*min(diff(t)),24*mean(diff(t)),24*max(diff(t)));
end
if max(t)-min(t)<365.25
    rd{3} = sprintf('    dur''n %#.3g dy; %s to %s GMT.',...
        max(t)-min(t),datestr(min(t),'dd-mmm-yy HH:MM'),...
        datestr(max(t),'dd-mmm-yy HH:MM'));
else
    rd{3} = sprintf('    dur''n %#.3g yr; %s to %s GMT.',...
        (max(t)-min(t))/365.25,datestr(min(t),'dd-mmm-yy HH:MM'),...
        datestr(max(t),'dd-mmm-yy HH:MM'));
end
if isequal(opt.method,'ols')
    rd{4} = '  Method: OLS.';
else
    rd{4} = sprintf(['  Method: IRLS, %s wt func; tunprm red''n fac '...
        '%#.3g.'],opt.method,opt.tunrdn);
end
if opt.linci
    rd{5} = '  Conf-ints: linearized;';
else
    rd{5} = sprintf('  Conf-ints: Monte-Carlo, %d realizations;',...
        opt.nrlzn);
end
if opt.white
    rd{5} = [rd{5} ' white.'];
else
    rd{5} = [rd{5} ' colored'];
    if ~opt.equi
        rd{5} = [rd{5} sprintf(' (Lmb-Scg ovrsmp %d).',...
            opt.lsfrqosmp)];
    else
        rd{5} = [rd{5} '.'];
    end
end
rd{6} = sprintf('  Model: allc=%d cnstits (%d non-ref), ',nNR+nR+nI,nNR);
if isequal(opt.cnstit,'auto')
    rd{6} = [rd{6} sprintf('auto (Rmin=%#.3g);',...
        opt.rmin)];
else
    rd{6} = [rd{6} 'specified as inputs;'];
end
if opt.notrend
    rd{7} = '    no trend included;';
else
    rd{7} = '    trend included;';
end
if isempty(opt.prefilt)
    rd{7} = [rd{7} ' no prefiltering correction;'];
else
    rd{7} = [rd{7} ' prefiltering correction applied;'];
end
if opt.nodsatnone
    rd{8} = '    no nod/sat corrections;';
elseif opt.nodsatlint
    rd{8} = sprintf(['    nod/sat corrections (lat=%.3f), '...
        'linearized times;'],lat);
else
    rd{8} = sprintf(['    exact nod/sat corrections (lat=%.3f);'],lat);
end
if opt.gwchnone
    rd{9} = '    g= raw phases (no Gwich correction);';
elseif opt.gwchlint
    rd{9} = '    g= Gwich phase lag, astr arg linearized times;';
else
    rd{9} = '    g= Gwich phase lag, astr arg exact times;';
end
if isempty(opt.infer)
    rd{10} = '    no inference.';
else
    rd{10} = '    exact inference:';
    if opt.inferaprx
        rd{10} = '    approximate inference:';
    end
    for j = 1:nI
        if opt.twodim
            rd{10+j} = sprintf(['     %s inferred from %s r+,r-= '...
                '%#.3g,%#.3g zta+,zta-= %#.3g,%#.3g deg.'],...
                deblank(opt.infer.infnam{j}),...
                deblank(opt.infer.refnam{j}),...
                opt.infer.amprat(j),opt.infer.amprat(nI+j),...
                opt.infer.phsoff(j),opt.infer.phsoff(nI+j));
        else
            rd{10+j} = sprintf(['     %s inferred from %s r= '...
                '%#.3g zta= %#.3g deg.'],...
                deblank(opt.infer.infnam{j}),...
                deblank(opt.infer.refnam{j}),...
                opt.infer.amprat(j),opt.infer.phsoff(j));
        end
    end
end
rundescr = char(rd);

%%--------------------------------------------------------- 
function E = ut_E(t,tref,frq,lind,lat,ngflgs,prefilt)
% UT_E()
% compute complex exponential basis function
% inputs
%   t = times [datenum UTC] (nt x 1)
%   tref = reference time [datenum UTC] (1 x 1)
%   frq = frequencies [cph] (nc x 1)
%   lind = list indices of constituents in ut_constants.mat (nc x 1)
%   lat = latitude [deg N] (1 x 1)
%   ngflgs = [NodsatLint NodsatNone GwchLint GwchNone] each 0/1
%       ([0 1 0 1] case not allowed, and not needed, in ut_E)
%   prefilt = 'prefilt' input to ut_solv
% output
%   E = complex exponential basis function [unitless] (nt x nc)
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 

nt = length(t);
nc = length(lind);
if ngflgs(2) && ngflgs(4)
    F = ones(nt,nc);
    U = zeros(nt,nc);
    V = 24*(t-tref)*frq';
else
    [F,U,V] = ut_FUV(t,tref,lind,lat,ngflgs);
end
E = F.*exp(1i*(U+V)*2*pi);
if ~isempty(prefilt)
    P=interp1(prefilt.frq,prefilt.P,frq)';
    P( P>max(prefilt.rng) | P<min(prefilt.rng) | isnan(P) )=1;
    E = E.*P(ones(nt,1),:);
end

%%--------------------------------------------------------- 
function [F,U,V] = ut_FUV(t,tref,lind,lat,ngflgs)
% UT_FUV()
% compute nodal/satellite correction factors and astronomical argument
% inputs
%   t = times [datenum UTC] (nt x 1)
%   tref = reference time [datenum UTC] (1 x 1)
%   lind = list indices of constituents in ut_constants.mat (nc x 1)
%   lat = latitude [deg N] (1 x 1)
%   ngflgs = [NodsatLint NodsatNone GwchLint GwchNone] each 0/1
% output
%   F = real nodsat correction to amplitude [unitless] (nt x nc)
%   U = nodsat correction to phase [cycles] (nt x nc)
%   V = astronomical argument [cycles] (nt x nc)
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 
% (uses parts of t_vuf.m from t_tide, Pawlowicz et al 2002)

nt = length(t);
nc = length(lind);
%% nodsat
if ngflgs(2) % none
    F = ones(nt,nc);
    U = zeros(nt,nc);
else
    if ngflgs(1) % linearized times
        tt = tref;
    else         % exact times
        tt = t;
    end
    ntt = length(tt);
    load('ut_constants.mat');
    [astro,~]=ut_astron(tt');
    if abs(lat)<5 
        lat=sign(lat).*5; 
    end
    slat=sin(pi*lat/180);
    rr=sat.amprat;
    j=find(sat.ilatfac==1);
    rr(j)=rr(j).*0.36309.*(1.0-5.0.*slat.*slat)./slat;
    j=find(sat.ilatfac==2);
    rr(j)=rr(j).*2.59808.*slat; 
    uu=rem( sat.deldood*astro(4:6,:)+sat.phcorr(:,ones(1,ntt)), 1);
    nfreq=length(const.isat); %#ok
    mat = rr(:,ones(1,ntt)).*exp(1i*2*pi*uu);
    F = ones(nfreq,ntt);
    ind = unique(sat.iconst);
    for i = 1:length(ind)
        F(ind(i),:) = 1+sum(mat(sat.iconst==ind(i),:),1);
    end
    U = imag(log(F))/(2*pi); % faster than angle(F)
    F=abs(F);
    for k=find(isfinite(const.ishallow))'
        ik=const.ishallow(k)+(0:const.nshallow(k)-1);
        j = shallow.iname(ik);
        exp1 = shallow.coef(ik);
        exp2 = abs(exp1);
        F(k,:)=prod(F(j,:).^exp2(:,ones(ntt,1)),1);
        U(k,:)=sum(U(j,:).*exp1(:,ones(ntt,1)),1);
    end
    F=F(lind,:)';
    U=U(lind,:)';
    if ngflgs(1) % nodal/satellite with linearized times
        F = F(ones(nt,1),:);
        U = U(ones(nt,1),:);
    end
end
%% gwch (astron arg)
if ngflgs(4) % none (raw phase lags not greenwich phase lags)
    if ~exist('const','var')
        load('ut_constants.mat','const');
    end
    [~,ader] = ut_astron(tref);
    ii=isfinite(const.ishallow); 
    const.freq(~ii) = (const.doodson(~ii,:)*ader)/(24);
    for k=find(ii)'
        ik=const.ishallow(k)+(0:const.nshallow(k)-1);
        const.freq(k)=sum(const.freq(shallow.iname(ik)).*shallow.coef(ik));
    end
    V = 24*(t-tref)*const.freq(lind)';
else 
    if ngflgs(3)  % linearized times
        tt = tref;
    else 
        tt = t;   % exact times
    end
    ntt = length(tt);
    if exist('astro','var')
        if ~isequal(size(astro,2),ntt)
            [astro,~]=ut_astron(tt');
        end        
    else
        [astro,~]=ut_astron(tt');
    end
    if ~exist('const','var')
        load('ut_constants.mat');
    end
    V=rem( const.doodson*astro+const.semi(:,ones(1,ntt)), 1);
    for k=find(isfinite(const.ishallow))'
        ik=const.ishallow(k)+(0:const.nshallow(k)-1);
        j = shallow.iname(ik);
        exp1 = shallow.coef(ik);
        V(k,:) = sum(V(j,:).*exp1(:,ones(ntt,1)),1);
    end
    V=V(lind,:)';
    if ngflgs(3)    % linearized times
        [~,ader] = ut_astron(tref);
        ii=isfinite(const.ishallow);
        const.freq(~ii) = (const.doodson(~ii,:)*ader)/(24);
        for k=find(ii)'
            ik=const.ishallow(k)+(0:const.nshallow(k)-1);
            const.freq(k)=sum( const.freq(shallow.iname(ik)).* ...
                shallow.coef(ik) );
        end
        V = V(ones(1,nt),:) + 24*(t-tref)*const.freq(lind)';
    end
end

%%--------------------------------------------------------- 
function [astro,ader] = ut_astron(jd)
% UT_ASTRON()
% calculate astronomical constants
% input
%   jd = time [datenum UTC] (1 x nt)
% outputs
%   astro = matrix [tau s h p np pp]T, units are [cycles] (6 x nt)
%   ader = matrix of derivatives of astro [cycles/day] (6 x nt)
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 
% (copy of t_astron.m from t_tide, Pawlowicz et al 2002)

d=jd(:)'-datenum(1899,12,31,12,0,0);
D=d/10000;
args=[ones(size(jd));
      d;
      D.*D;
      D.^3];
sc= [ 270.434164,13.1763965268,-0.0000850, 0.000000039];
hc= [ 279.696678, 0.9856473354, 0.00002267,0.000000000];
pc= [ 334.329556, 0.1114040803,-0.0007739,-0.00000026];
npc=[-259.183275, 0.0529539222,-0.0001557,-0.000000050];
ppc=[ 281.220844, 0.0000470684, 0.0000339, 0.000000070];
astro=rem( [sc;hc;pc;npc;ppc]*args./360.0 ,1);
tau=rem(jd(:)',1)+astro(2,:)-astro(1,:);
astro=[tau;astro];
dargs=[zeros(size(jd));
       ones(size(jd));
       2.0e-4.*D;
       3.0e-4.*D.*D];
ader=[sc;hc;pc;npc;ppc]*dargs./360.0;
dtau=1.0+ader(2,:)-ader(1,:);
ader=[dtau;ader];

%%--------------------------------------------------------- 
function [Lsmaj,Lsmin,theta,g] = ut_cs2cep(XY)
% UT_CS2CEP()
% compute current ellipse parameters from cosine-sine coefficients
% inputs 
%   two-dim case: XY = [Xu Yu Xv Yv] 4-column matrix 
%   one-dim case: XY = [Xu Yu] 2-column matrix 
%                      OR 4-column w/ Xv=Yv=zeros(size(Xu))
%      where: Xu,Yu are cosine, sine coeffs of u, & same for v
% outputs
%   two-dim case:
%     Lsmaj, Lsmin = column vectors [units of XY] (size of Xu)
%     theta = column vector [deg. ccw rel. +x-axis, 0-180] (size of Xu)
%     g = column vector [degrees, 0-360] (size of Xu)
%   one-dim case: same, where Lsmaj = A, and Lsmin and theta are zeros
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 

Xu = XY(:,1);
Yu = XY(:,2);
if size(XY,2)>2
    Xv = XY(:,3);
    Yv = XY(:,4);
else
    Xv = zeros(size(Xu));
    Yv = zeros(size(Yu));
end
ap = ((Xu+Yv)+1i*(Xv-Yu))/2;
am = ((Xu-Yv)+1i*(Xv+Yu))/2;
Ap = abs(ap);
Am = abs(am);
Lsmaj = Ap+Am;
Lsmin = Ap-Am;
epsp = angle(ap)*180/pi;
epsm = angle(am)*180/pi;
theta = mod((epsp+epsm)/2,180);
g = mod(-epsp+theta,360);

%%--------------------------------------------------------- 
function P = ut_pdgm(t,e,cfrq,equi,frqosmp)
% UT_PDGM()
% spectral densities (band-averages of line-decimated periodograms)
% inputs
%   t = times [datenum UTC] (nt x 1)
%   e = error (residual; model-raw misfit), complex/real 2D/1D (nt x 1)
%   cfrq = frequencies of NR & R constituents [cycles per hour] (nc x 1)
%   equi = 0/1 if equispaced times or not (1 x 1)
%   frqosamp = lomb-scargle freq oversampling factor (ignored if equi=1)
% outputs (all 9 x 1)
%   P.Puu= 1-sided auto-sp dens of u err (=real(e)) [e units^2/cph] (9 x 1)
%   P.fbnd= edges of freq bands [cph] (9 x 2) 
%   P.Pvv (2dim case only)= as P.Puu but for v err (=imag(e)) (9 x 1) 
%   P.Puv (2dim case only)= cross-sp dens between u and v (9 x 1)
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 
% (band-averaging and line decimating from t_tide, Pawlowicz et al 2002)

nt = length(e);
hn = hanning(nt);
if equi
    [Puu1s,allfrq] = pwelch(real(e),hn,0,nt);
else
    [Puu1s,allfrq] = ut_lmbscga(real(e),t,hn,frqosmp);
end
fac = (nt-1)/(2*pi*(t(end)-t(1))*24); % conv fac: rad/sample to cph
allfrq = allfrq*fac; % to [cycle/hour] from [rad/samp]
Puu1s = Puu1s/fac; % to [e units^2/cph] from [e units^2/(rad/samp)]
[P.Puu,P.fbnd] = ut_fbndavg(Puu1s,allfrq,cfrq);    
if ~isreal(e)
    if equi
        [Pvv1s,~] = pwelch(imag(e),hn,0,nt);
        [Puv1s,~] = cpsd(real(e),imag(e),hn,0,nt);
    else
        [Pvv1s,~] = ut_lmbscga(imag(e),t,hn,frqosmp);
        [Puv1s,~] = ut_lmbscgc(real(e),imag(e),t,hn,frqosmp);
    end
    Pvv1s = Pvv1s/fac;
    [P.Pvv,~] = ut_fbndavg(Pvv1s,allfrq,cfrq);
    Puv1s = real(Puv1s)/fac;
    [P.Puv,~] = ut_fbndavg(Puv1s,allfrq,cfrq);
    P.Puv = abs(P.Puv);
end

%%--------------------------------------------------------- 
function [Pxx,F] = ut_lmbscga(x,t,w,ofac)
% UT_LMBSCA()
% Auto-spectral density estimate based on the mean-removed Lomb-Scargle 
%   periodogram; designed such that, in the special case of equispaced
%   input times and no frequency oversampling, the output units,
%   normalization, and frequency values replicate those of pwelch(). 
%
% Inputs: 
%   * Sequence of real values x from n arbitrarily distributed times t
%   * Taper/window shape w at n equispaced times (e.g. w=hanning(n) )
%   * Frequency oversampling factor ofac (integer >=1)
%
% Outputs:
%   * One-sided auto-spectral density estimate Pxx [units: (x-units^2) /
%      (radians per sample)] as column vector, based on mean-removed and
%      unnormalized Lomb-Scargle periodogram, after applying window/taper w
%      to entire record (segmenting is not implemented); at zero and
%      Nyquist frequencies Pxx is computed by the FFT-based auto-spectral
%      density formula.
%   * Frequencies F (units: radians per sample) of the Pxx estimates, as 
%      column vector: for ofac=1, these are the same equispaced values as
%      for an FFT-based spectral calculation using a time series of n
%      values at uniformly distributed times from min(t) to max(t), thus F
%      has resolution 2*pi/n (rad/sample) and for even n it spans the
%      interval [0,Pi] with (n+2)/2 values, while for odd n it spans the
%      interval [0,Pi) with (n+1)/2 values; for ofac>1, the same ranges are
%      spanned but with resolution increased by ofac times.
%
% Notes:
%   Designed such that, for the special case of input values x from 
%      uniformly distributed times t, [Pxx,F]=ut_lmbscga(x,t,w,1) 
%      gives results (both Pxx and F) identical to [Pxx,F]=pwelch(x,w,0,n).
%   The cross-periodogram companion function is ut_lmbscgc().
%   Finally, for context with respect to the now-deprecated psd() 
%      function, the functions pwelch() and psd() are related as follows.
%      Consider [Pxxw,Fw] = pwelch(x,w,0,n) and [Pxxp,Fp] = psd(x,n,1/dt),
%      for real x, w = hanning(n), and sampling interval dt [t units]. 
%      A. Pxxw is the one-sided power spectral density with units 
%      [x-units^2/(rad/sample)] and Pxxp is the two-sided power spectrum
%      (not spectral density) with units [x-units^2]. Therefore, to obtain
%      Pxxw from Pxxp requires (i) division by the sampling frequency
%      (1/dt) because Pxxp is power spectrum not power spectral density,
%      yielding dt*Pxxp [units: x-units^2/(cycles per t-unit)]; (ii)
%      division by 2*pi and by dt, to convert from cycles per t-unit to
%      radians per sample, and (iii) multiplication by 2 to convert from
%      two-sided to one-sided, for frequencies other than 0 and Nyquist. 
%       Even n: Pxxw is [Pxxp(1)/(2*pi) Pxxp(2:end-1)/pi Pxxp(end)/(2*pi)];
%       Odd n:  Pxxw is [Pxxp(1)/(2*pi) Pxxp(2:end)/pi].
%      B. Units of Fw are radians per sample and units of Fp are cycles 
%      per t-unit, so to obtain Fw from Fp requires multiplying by 2*pi 
%      and by dt. Even or odd n: Fw is 2*pi*dt*Fp.
% UTide v1p0 9/2011 d.codiga@gso.uri.edu  

% reshape and check inputs
x = x(:)'; 
t = t(:)';
w = w(:)';
if ~isequal(size(x),size(t),size(w)) || ~(isreal(x)&&isreal(t)&&isreal(w)) 
    error('ut_lmbscga: x,t, and w must be same size and real.');
end
ofac = round(ofac);
if ofac<1
    error('ut_lmbscga: ofac must be >= 1.');
end
n = length(x);
% taper
dt = (max(t) - min(t))/(n-1);
w = interp1(min(t):dt:max(t),w,t);
xw = x.*w;
% frequencies (excluding 0 and nyquist)
F = ((1/ofac):(1/ofac):n/2)'*2*pi/n; % radian/sample
if ~mod(n,2) % even, exclude nyqyuist
    F(end) = [];
end
% unnormalized, 2-sided, mean-removed lomb-scargle periodogram [x-units^2] 
arg = (2/dt)*(F*t)';
tau = 0.5*imag(log(sum(cos(arg))+1i*sum(sin(arg)))); % see angle() hdr
arg = 0.5*arg - tau(ones(n,1),:);
xb = mean(xw);
Pxx = sum(sparse(1:n,1:n,xw-xb)*cos(arg)).^2./sum(cos(arg).^2);
Pxx = 0.5*(Pxx+sum(sparse(1:n,1:n,xw-xb)*sin(arg)).^2./sum(sin(arg).^2));
% include zero-frequency
Pxx = [n*xb^2; Pxx']; 
F = [0; F];
if ~mod(n,2)  % even, include nyquist
    c = [ones(1,n/2); -1*ones(1,n/2)];
    Pxx = [Pxx; n*mean(xw.*c(:)')^2];
    F = [F; pi];
end
% recover taper suppression
Wss = n*sum(w.*w);
Pxx = n^2*Pxx/Wss;
% divide by 2*pi to get spectral density rel radians/sample
Pxx = Pxx/(2*pi); % two-sided, [x-units^2/(rad/sample)]
% multiply by 2 to get one-sided (except zero and nyquist)
Pxx(2:end-1) = 2*Pxx(2:end-1); % one-sided, [x-units^2/(rad/sample)]
if mod(n,2) % odd n, highest frequency not nyquist
    Pxx(end) = 2*Pxx(end); % one-sided, [x-units^2/(rad/sample)]
end

%%--------------------------------------------------------- 
function [Pxy,F] = ut_lmbscgc(x,y,t,w,ofac)
% UT_LMBSCGC()
% Cross-spectral density estimate based on the mean-removed Lomb-Scargle 
%   cross-periodogram; designed such that, in special case of equispaced
%   input times and no frequency oversampling, the output units,
%   normalization, and frequency values replicate those of CPSD().  
%
% Inputs: 
%   * Real-valued sequences x and y from n arbitrarily distributed times t
%   * Taper/window shape w at n equispaced times (e.g. w=hanning(n) )
%   * Frequency oversampling factor ofac (integer >=1)
%
% Outputs:
%   * One-sided cross-spectral density estimate Pxy [units: (x-units)
%      (y-units)/(radians per sample)], based on the mean-removed and
%      unnormalized Lomb-Scargle cross-periodogram, after applying
%      window/taper w to the entire records (segmenting is not
%      implemented); at zero and Nyquist frequencies Pxy is computed by the
%      FFT-based cross-spectral density formula.
%   * Frequencies F (units: radians per sample) of the Pxy estimates: for 
%      ofac=1, these are the same equispaced values as for an FFT-based 
%      spectral calculation using a time series of n values at uniformly
%      distributed times from min(t) to max(t), thus F has resolution
%      2*pi/n (rad/sample) and for even n it spans the interval [0,Pi] with
%      (n+2)/2 values, while for odd n it spans the interval [0,Pi) with
%      (n+1)/2 values; for ofac>1, the same ranges are spanned but with
%      resolution increased by ofac times.
%
% Notes:
%   Designed such that, for the special case of input values x and y from 
%      uniformly distributed times t, [Pxy,F]=ut_lmbscgc(x,y,t,w,1) 
%      gives results (both Pxy and F) identical to [Pxx,F]=cpsd(x,y,w,0,n).
%   The auto-periodogram companion function is ut_lmbscga(). 
%   The relationship between cpsd() and the deprecated csd() function is
%      fully analogous to the relationship between pwelch() and the
%      deprecated psd() function; an explanation of the latter is given at
%      the end of the function description text in ut_lmbscga.m.
% UTide v1p0 9/2011 d.codiga@gso.uri.edu

% reshape and check inputs
x = x(:)';
y = y(:)';
t = t(:)';
w = w(:)';
if ~isequal(size(x),size(y),size(t),size(w)) || ...
        ~(isreal(x)&&isreal(y)&&isreal(t)&&isreal(w)) 
    error('ut_lmbscgc: x,y,t, and w must be same size and real.');
end
ofac = round(ofac);
if ofac<1
    error('ut_lmbscgc: ofac must be >= 1.');
end
n = length(x);
% taper
dt = (max(t) - min(t))/(n-1);
w = interp1(min(t):dt:max(t),w,t);
xw = x.*w;
yw = y.*w;
% frequencies (excluding 0 and nyquist)
F = ((1/ofac):(1/ofac):n/2)'*2*pi/n; % radian/sample
if ~mod(n,2) % even, exclude nyqyuist
    F(end) = [];
end
% unnormalized, 2-sided, mean-removed lomb-scargle cross-periodogram 
% [(x-units)(y-units)] 
arg = (1/dt)*(F*t)';
tau = 0.5*imag(log( sum(cos(2*arg)) + 1i*sum(sin(2*arg)) ));
xb = mean(xw);
yb = mean(yw);
X1 = sum(sparse(1:n,1:n,xw-xb)*cos(arg));
X2 = sum(sparse(1:n,1:n,xw-xb)*sin(arg));
Y1 = sum(sparse(1:n,1:n,yw-yb)*cos(arg));
Y2 = sum(sparse(1:n,1:n,yw-yb)*sin(arg));
arg = arg - tau(ones(n,1),:);
A = 1./sqrt(sum(cos(arg).^2));
B = 1./sqrt(sum(sin(arg).^2));
F0=exp(1i*tau);
X = F0.*(A.*X1 + 1i*B.*X2);
Y = conj(F0).*(A.*Y1 - 1i*B.*Y2);
Pxy = 0.5*X.*Y;
% include zero-frequency
Pxy = [n*xb*yb; Pxy'];
F = [0; F];
if ~mod(n,2)  % even, include nyquist freq
    c = [ones(1,n/2); -1*ones(1,n/2)];
    Pxy = [Pxy; n*mean(xw.*c(:)')*mean(yw.*c(:)')];
    F = [F; pi];
end
% recover taper suppression
Wss = n*sum(w.*w);
Pxy = n^2*Pxy/Wss;
% divide by 2*pi to get spectral density rel radians/sample
Pxy = Pxy/(2*pi); % two-sided, [x-units^2/(rad/sample)]
% multiply by 2 to get one-sided (except zero and nyquist)
Pxy(2:end-1) = 2*Pxy(2:end-1); % one-sided, [x-units^2/(rad/sample)]
if mod(n,2) % odd n, highest frequency not nyquist
    Pxy(end) = 2*Pxy(end); % one-sided, [x-units^2/(rad/sample)]
end

%%--------------------------------------------------------- 
function [avP,fbnd] = ut_fbndavg(P,allfrq,cfrq)
% UT_FBNDAVG()
% line-decimate and band-average spectra
% inputs
%   P = periodogram to treat [e units^2/cph]
%   allfrq = frequency values of (equispaced) P estimates [cph]
%   cfrq = frequencies of constituents [cph] (nc x 1)
% outputs
%   avP = line-decimated and band-averaged spectrum [e units^2/cph] (9 x 1)
%   fbnd = frequencies [cph] at edges of averaged bands (9 x 2)
% UTide v1p0 9/2011 d.codiga@gso.uri.edu
% (based on residual_spectrum.m of t_tide, Pawlowicz et al 2002)

df=allfrq(3)-allfrq(2);
P(round(cfrq./df)+1)=NaN; 
fbnd =[.00010 .00417;
    .03192 .04859;
    .07218 .08884;
    .11243 .12910;
    .15269 .16936;
    .19295 .20961;
    .23320 .25100;
    .26000 .29000;
    .30000 .50000];
nfbnd=size(fbnd,1);
avP=zeros(nfbnd,1);
for k=nfbnd:-1:1,
    b1 = find(allfrq>=fbnd(k,1));
    b2 = find(allfrq<=fbnd(k,2));
    b3 = find(isfinite(P));
    jbnd=intersect(intersect(b1,b2),b3); 
    if any(jbnd),
        avP(k)=mean(P(jbnd));
    elseif k<nfbnd,
        avP(k)=P(k+1);   
    end
end

%%--------------------------------------------------------- 
function [sig1,sig2] = ut_linci(X,Y,sigX,sigY)
% UT_LINCI()
% current ellipse parameter uncertainties from cosine/sine coefficient
%   uncertainties, by linearized relations w/ correlations presumed zero
% inputs: (two-dim case complex, one-dim case real)
%   X = Xu + i*Xv 
%   Y = Yu + i*Yv
%       for Xu =real(X) = u cosine coeff; Yu =real(Y) = u sine coeff
%           Xv =imag(X) = v cosine coeff; Yv =imag(Y) = v sine coeff
%   sigX = sigXu + i*sigXv
%   sigY = sigYu + i*sigYv
%       for sigXu =real(sigX) =stddev(Xu); sigYu =real(sigY) =stddev(Yu)
%           sigXv =imag(sigX) =stddev(Xv); sigYv =imag(sigY) =stddev(Yv)
% outputs: 
%   two-dim case, complex
%       sig1 = sig_Lsmaj +1i*sig_Lsmin [same units as inputs]
%       sig2 = sig_g + 1i*sig_theta [degrees]
%   one-dim case, real 
%       sig1 = sig_A [same units as inputs]
%       sig2 = sig_g [degrees]
% UTide v1p0 9/2011 d.codiga@gso.uri.edu
% (adapted from errell.m of t_tide, Pawlowicz et al 2002)

Xu = real(X(:));
sigXu = real(sigX);
Yu = real(Y(:));
sigYu = real(sigY);

Xv = imag(X(:));
sigXv = imag(sigX(:));
Yv = imag(Y(:));
sigYv = imag(sigY(:));

rp=.5.*sqrt((Xu+Yv).^2+(Xv-Yu).^2);
rm=.5.*sqrt((Xu-Yv).^2+(Xv+Yu).^2);
sigXu2=sigXu.^2;sigYu2=sigYu.^2;
sigXv2=sigXv.^2;sigYv2=sigYv.^2;

ex=(Xu+Yv)./rp;
fx=(Xu-Yv)./rm;
gx=(Yu-Xv)./rp;
hx=(Yu+Xv)./rm;

% major axis
dXu2=(.25.*(ex+fx)).^2;
dYu2=(.25.*(gx+hx)).^2;
dXv2=(.25.*(hx-gx)).^2;
dYv2=(.25.*(ex-fx)).^2;
sig1 = sqrt(dXu2.*sigXu2+dYu2.*sigYu2+dXv2.*sigXv2+dYv2.*sigYv2);

% phase 
rn=2.*(Xu.*Yu+Xv.*Yv);
rd=Xu.^2-Yu.^2+Xv.^2-Yv.^2;
den=rn.^2+rd.^2;
dXu2=((rd.*Yu-rn.*Xu)./den).^2;
dYu2=((rd.*Xu+rn.*Yu)./den).^2;
dXv2=((rd.*Yv-rn.*Xv)./den).^2;
dYv2=((rd.*Xv+rn.*Yv)./den).^2;
sig2 = (180/pi)*sqrt(dXu2.*sigXu2+dYu2.*sigYu2+dXv2.*sigXv2+dYv2.*sigYv2);

if ~isreal(X)
    % minor axis 
    dXu2=(.25.*(ex-fx)).^2;
    dYu2=(.25.*(gx-hx)).^2;
    dXv2=(.25.*(hx+gx)).^2;
    dYv2=(.25.*(ex+fx)).^2;
    sig1 = sig1 + ...
        1i*sqrt(dXu2.*sigXu2+dYu2.*sigYu2+dXv2.*sigXv2+dYv2.*sigYv2);
    
    % orientation
    rn=2.*(Xu.*Xv+Yu.*Yv);
    rd=Xu.^2+Yu.^2-(Xv.^2+Yv.^2);
    den=rn.^2+rd.^2;
    dXu2=((rd.*Xv-rn.*Xu)./den).^2;
    dYu2=((rd.*Yv-rn.*Yu)./den).^2;
    dXv2=((rd.*Xu+rn.*Xv)./den).^2;
    dYv2=((rd.*Yu+rn.*Yv)./den).^2;
    sig2 = sig2 + 1i*(180/pi).*sqrt(dXu2.*sigXu2+dYu2.*sigYu2 + ...
        dXv2.*sigXv2+dYv2.*sigYv2);    
end

%%--------------------------------------------------------- 
function covposdef = ut_nearposdef(cov,maxit)
% UT_NEARPOSDEF()
% find nearest positive definite matrix by Higham 2002 algorithm
%   with identity weight matrix
% inputs
%   cov - matrix for which to find the nearest positive definite matrix
%   maxit - max iterations before quitting (optional; default 1000)
% output
%   covposdef - the nearest positive definite covariance matrix
%             - diag(diag(cov)) if max iterations limit is met 
% UTide v1p0 9/2011 d.codiga@gso.uri.edu
% (builds on validcorr.m by Erland Ringstad, by adding:
%   conversion from correlation matrix to covariance matrix;
%   end of search once all eigenvalues are positive;
%   optional maxit input. )

if size(cov,1)~=size(cov,2)
    fprintf('ut_nearposdef: cov must be square and real.\n');
    covposdef = nan;
    return;
end
if all(eig(cov)>=0)
    covposdef = cov;
    return;
end
if ~exist('maxit','var')
    maxit = 1000;
elseif maxit<=0 || rem(maxit,1)
    fprintf('ut_nearposdef: maxit must be positive integer.\n');
    covposdef = nan;
    return;
end
S = zeros(size(cov));
Y = cov;
k = 0;
while  ~all(eig(Y)>=0) && k<maxit
    R = Y - S;
    [Q,D] = eig(R);
    X = Q*max(D,0)*Q';
    S = X - R;
    Y = X - diag(diag(X)-1);
    k = k + 1;
end
if ~all(eig(Y)>=0)
    fprintf(sprintf(['ut_nearposdef: reached max iterations'...
        ' (%d); zeroing off-diagonals.\n'],maxit));
    covposdef = diag(diag(cov));
else
    std = sqrt(diag(cov));
    covposdef = Y.*(std*std');
end

%%--------------------------------------------------------- 
function ain = ut_cluster(ain,clusang)
% UT_CLUSTER()
% UTide v1p0 9/2011 d.codiga@gso.uri.edu
% (copy of t_cluster.m of t_tide, Pawlowicz et al 2002)
ii=(ain-ain(:,ones(1,size(ain,2))))>clusang/2;
ain(ii)=ain(ii)-clusang;
ii=(ain-ain(:,ones(1,size(ain,2))))<-clusang/2;
ain(ii)=ain(ii)+clusang;

%%--------------------------------------------------------- 
function [diagn,usnrc,vsnrc] = ut_diagntable(coef,cnstit,t,u,v,...
    xmod,m,B,W,varMSM,Gall,Hall,elor,varcov_mC,indPE)
% UT_DIAGNTABLE()
% compute and store diagnostics, then create character array table
% inputs: see ut_solv1()
% outputs:
%   diagn = diagnostic quantities results (see report)
%   usnrc, vsnrc = reconstructed superposed harmonics using constituents
%       with SNR above the threshold
% UTide v1p0 9/2011 d.codiga@gso.uri.edu

nNR = length(cnstit.NR.frq);
nR = 0;
nI = 0;
if ~isempty(coef.aux.opt.infer)
    nR = length(cnstit.R);
    nI = length(coef.aux.opt.infer.infnam);
end    
nallc = nNR+nR+nI;
nc = nNR+nR;
diagn = coef.diagn;
tref = coef.aux.reftime;
[usnrc,vsnrc] = ut_diagnrcn(t,coef,diagn.SNR,...
    coef.aux.opt.diagnminsnr);
% TVallc, TVraw, PTVallc
if coef.aux.opt.twodim
    urawtid = u - coef.umean;
    vrawtid = v - coef.vmean;
    umodtid = real(xmod) - coef.umean;
    vmodtid = imag(xmod) - coef.vmean;
    usnrctid = usnrc - coef.umean;
    vsnrctid = vsnrc - coef.vmean;
    if ~coef.aux.opt.notrend
        umodtid = umodtid - coef.uslope*(t-tref);
        vmodtid = vmodtid - coef.vslope*(t-tref);
        urawtid = urawtid - coef.uslope*(t-tref);
        vrawtid = vrawtid - coef.vslope*(t-tref);
        usnrctid = usnrctid - coef.uslope*(t-tref);
        vsnrctid = vsnrctid - coef.vslope*(t-tref);
    end
    diagn.TVraw = mean(urawtid.*urawtid + vrawtid.*vrawtid);
    diagn.TVallc = mean(umodtid.*umodtid + vmodtid.*vmodtid);
    diagn.TVsnrc = mean(usnrctid.*usnrctid + vsnrctid.*vsnrctid);
else
    urawtid = u - coef.mean;
    umodtid = real(xmod) - coef.mean;
    usnrctid = usnrc - coef.mean;
    if ~coef.aux.opt.notrend
        urawtid = urawtid - coef.slope*(t-tref);
        umodtid = umodtid - coef.slope*(t-tref);
        usnrctid = usnrctid - coef.slope*(t-tref);
    end
    diagn.TVraw = mean(urawtid.*urawtid);
    diagn.TVallc = mean(umodtid.*umodtid);
    diagn.TVsnrc = mean(usnrctid.*usnrctid);
end
diagn.PTVallc = 100*diagn.TVallc/diagn.TVraw;
diagn.PTVsnrc = 100*diagn.TVsnrc/diagn.TVraw;
% SNRallc, K
diagn.SNRallc = real(((m'*B')*W*(B*m))/varMSM);
diagn.K = cond(B);
% RR, RNM, CorMx
diagn.lo.name = cellstr(repmat('none',nallc,1));
diagn.lo.RR = nan*ones(nallc,1);
diagn.lo.RNM = diagn.lo.RR;
diagn.lo.CorMx = diagn.lo.RR;
diagn.hi = diagn.lo;
[~,ind] = sort(coef.aux.frq(1:nc),'ascend');
for i = 1:nc-1
    c1 = ind(i);
    c2 = ind(i+1);
    diagn.lo.name(i+1) = coef.name(c1);
    diagn.hi.name(i) = coef.name(c2);
    diagn.lo.RR(i+1) = 24*elor*(coef.aux.frq(c2) - ...
        coef.aux.frq(c1))/coef.aux.opt.rmin;
    diagn.lo.RNM(i+1) = diagn.lo.RR(i+1)* ...
        sqrt((diagn.SNR(c1)+diagn.SNR(c2))/2);
    den = sqrt(diag(squeeze(varcov_mC(c1,:,:)))* ...
        diag(squeeze(varcov_mC(c2,:,:)))');
    G =[Gall(c1,c2) Gall(c1,c2+nc); Gall(c1+nc,c2) Gall(c1+nc,c2+nc);];
    H =[Hall(c1,c2) Hall(c1,c2+nc); Hall(c1+nc,c2) Hall(c1+nc,c2+nc);];
    cDuu = nan*ones(2,2);
    if ~coef.aux.opt.twodim
        cDuu(1,1) = real(sum(sum(G)))/2;
        cDuu(1,2) = imag(H(1,1)-H(1,2)+H(2,1)-H(2,2))/2;
        cDuu(2,1) = imag(-G(1,1)-G(1,2)+G(2,1)+G(2,2))/2;
        cDuu(2,2) = real(H(1,1)-H(1,2)-H(2,1)+H(2,2))/2;
        diagn.lo.CorMx(i+1) = max(max(abs(cDuu)./den));
    else
        cDuv(1,1) = imag(sum(sum(-H)))/2;
        cDuv(1,2) = real(G(1,1)-G(1,2)+G(2,1)-G(2,2))/2;
        cDuv(2,1) = real(-H(1,1)-H(1,2)+H(2,1)+H(2,2))/2;
        cDuv(2,2) = imag(-G(1,1)+G(1,2)+G(2,1)-G(2,2))/2;
        cDvu(1,1) = imag(sum(sum(G)))/2;
        cDvu(1,2) = real(-H(1,1)+H(1,2)-H(2,1)+H(2,2))/2;
        cDvu(2,1) = real(-G(1,1)-G(1,2)+G(2,1)+G(2,2))/2;
        cDvu(2,2) = imag(H(1,1)-H(1,2)-H(2,1)+H(2,2))/2;
        cDvv(1,1) = real(sum(sum(H)))/2;
        cDvv(1,2) = imag(G(1,1)-G(1,2)+G(2,1)-G(2,2))/2;
        cDvv(2,1) = imag(-H(1,1)-H(1,2)+H(2,1)+H(2,2))/2;
        cDvv(2,2) = real(G(1,1)-G(1,2)-G(2,1)+G(2,2))/2;
        cv = [cDuu cDuv; cDvu cDvv;];
        diagn.lo.CorMx(i+1) = max(max(abs(cv)./den));
    end
end
diagn.hi.RR(1:nc-1) = diagn.lo.RR(2:nc);
diagn.hi.RNM(1:nc-1) = diagn.lo.RNM(2:nc);
diagn.hi.CorMx(1:nc-1) = diagn.lo.CorMx(2:nc);
% sort
ind(ind) = 1:length(ind);
diagn.lo.name(1:nc) = diagn.lo.name(ind);
diagn.lo.RR(1:nc) = diagn.lo.RR(ind);
diagn.lo.RNM(1:nc) = diagn.lo.RNM(ind);
diagn.lo.CorMx(1:nc) = diagn.lo.CorMx(ind);
diagn.hi.name(1:nc) = diagn.hi.name(ind);
diagn.hi.RR(1:nc) = diagn.hi.RR(ind);
diagn.hi.RNM(1:nc) = diagn.hi.RNM(ind);
diagn.hi.CorMx(1:nc) = diagn.hi.CorMx(ind);
diagn.lo.name = diagn.lo.name(indPE);
diagn.lo.RR = diagn.lo.RR(indPE);
diagn.lo.RNM = diagn.lo.RNM(indPE);
diagn.lo.CorMx = diagn.lo.CorMx(indPE);
diagn.hi.name = diagn.hi.name(indPE);
diagn.hi.RR = diagn.hi.RR(indPE);
diagn.hi.RNM = diagn.hi.RNM(indPE);
diagn.hi.CorMx = diagn.hi.CorMx(indPE);
diagn.SNR = diagn.SNR(indPE);

% summary table char array
minsnr = coef.aux.opt.diagnminsnr;
tbl = cell(nallc+6+sign(nI)*(nI+1),1);
tbl{1} = 'UTide Summary Diagnostics Table:';
tbl{2} = sprintf(['   Rmin= %#-9.3g   MinSNR= %#-9.3g   '...
    '(* SNR >= MinSNR)'],coef.aux.opt.rmin,minsnr);
tbl{3} = sprintf('   K= %#-9.3g   SNRallc= %#-9.3g',...
    diagn.K,diagn.SNRallc);
tbl{4} = sprintf(['   TVallc= %#-9.3g    TVsnrc= %#-9.3g   '...
    'TVraw= %#-9.3g'],diagn.TVallc,diagn.TVsnrc,...
    diagn.TVraw);
tbl{5} = sprintf('   PTVallc= %5.1f%%   PTVsnrc= %5.1f%%',...
    diagn.PTVallc,diagn.PTVsnrc);
tbl{6} = [' NAME    PE       SNR loNAME   loRR    loRNM  '...
    'loCorMx hiNAME   hiRR    hiRNM  hiCorMx'];
main = [blanks(nallc)' char(diagn.name)];
main(diagn.SNR>=minsnr,1) = repmat('*',sum(diagn.SNR>=minsnr),1);
main = [main reshape(sprintf('%6.2f%%',diagn.PE),7,nallc)'];
main = [main reshape(sprintf(' %#8.2g',diagn.SNR),9,nallc)'];
if nI
    indPE(indPE) = 1:nallc;
    diagn.lo.name(indPE(end-nI+1:end)) = ...
        cellstr(repmat('(I) ',nI,1));
    diagn.hi.name(indPE(end-nI+1:end)) = ...
        cellstr(repmat('(I) ',nI,1));
end
lo = char(diagn.lo.name);
lo = [lo reshape(sprintf(' %#8.2g',diagn.lo.RR),9,nallc)'];
lo = [lo reshape(sprintf(' %#8.2g',diagn.lo.RNM),9,nallc)'];
lo = [lo reshape(sprintf(' %#8.2g',diagn.lo.CorMx),9,nallc)'];
hi = char(diagn.hi.name);
hi = [hi reshape(sprintf(' %#8.2g',diagn.hi.RR),9,nallc)'];
hi = [hi reshape(sprintf(' %#8.2g',diagn.hi.RNM),9,nallc)'];
hi = [hi reshape(sprintf(' %#8.2g',diagn.hi.CorMx),9,nallc)'];
tbl(7:7+nallc-1) = cellstr([main repmat(' ',nallc,1) lo ...
    repmat(' ',nallc,1) hi]);
if nI
    tbl{end-nI} = '  Inferred constituents (I):';
    cntr =0;
    for i=1:nR
        for j = 1:cnstit.R{i}.nI
            cntr = cntr+1;
            if coef.aux.opt.twodim
                tbl{end-nI+cntr} = sprintf(['    %s inferred from'...
                    ' reference %s (r+,r- = %#.3g,%#.3g; ',...
                    'zta+,zta- = %#.3g,%#.3g deg)'],...
                    deblank(cnstit.R{i}.I.name{j}),...
                    deblank(cnstit.R{i}.name),...
                    abs(cnstit.R{i}.I.Rp(j)),...
                    abs(cnstit.R{i}.I.Rm(j)),...
                    (180/pi)*angle(cnstit.R{i}.I.Rp(j)),...
                    -(180/pi)*angle(cnstit.R{i}.I.Rm(j)));
            else
                tbl{end-nI+cntr} = sprintf(['    %s inferred from'...
                    ' reference %s (r = %#.3g zta = %#.3g deg)'],...
                    deblank(cnstit.R{i}.I.name{j}),...
                    deblank(cnstit.R{i}.name),...
                    abs(cnstit.R{i}.I.Rp(j)),...
                    (180/pi)*angle(cnstit.R{i}.I.Rp(j)));
            end
        end
    end
end
diagn.table = char(reshape(strrep(tbl(:)','NaN','---')',...
    size(tbl,1),size(tbl,2)));

%%--------------------------------------------------------- 
function [u,v] = ut_diagnrcn(t,coef,SNR,MinSNR)
% UT_DIAGNRCN()
% abbreviated version of ut_reconstr1 solely for use by ut_diagntable:
%       accepts only non-nan times
%       instead of varargin, requires that you pass in
%           SNR (in NR-R-I order, same as coef.Lsmaj etc, or coef.A etc) 
%           and MinSNR
% does not allow for subsetting by PE or specified constits
%   (as ut_reconstr1 does)
% see comments for ut_reconstr()
% UTide v1p0 9/2011 d.codiga@gso.uri.edu

ind = SNR>=MinSNR;
rpd = pi/180;
% complex coefficients
if coef.aux.opt.twodim
    ap = 0.5*(coef.Lsmaj(ind) + coef.Lsmin(ind)) .* ...
        exp(1i*(coef.theta(ind) - coef.g(ind))*rpd);
    am = 0.5*(coef.Lsmaj(ind) - coef.Lsmin(ind)) .* ...
        exp(1i*(coef.theta(ind) + coef.g(ind))*rpd);
else
    ap = 0.5*coef.A(ind).*exp(-1i*coef.g(ind)*rpd);
    am = conj(ap);
end
% exponentials
ngflgs = [coef.aux.opt.nodsatlint coef.aux.opt.nodsatnone ...
    coef.aux.opt.gwchlint coef.aux.opt.gwchnone];
E = ut_E(t,coef.aux.reftime,coef.aux.frq(ind),coef.aux.lind(ind),...
    coef.aux.lat,ngflgs,coef.aux.opt.prefilt);
% fit calc
fit = E*ap + conj(E)*am;
% mean (& trend)
if coef.aux.opt.twodim
    if coef.aux.opt.notrend
        u = real(fit) + coef.umean;
        v = imag(fit) + coef.vmean;
    else
        u = real(fit) + coef.umean + ...
            coef.uslope*(t-coef.aux.reftime);
        v = imag(fit) + coef.vmean + ...
            coef.vslope*(t-coef.aux.reftime);
    end
else
    if coef.aux.opt.notrend
        u = real(fit) + coef.mean;
    else
        u = real(fit) + coef.mean + ...
            coef.slope*(t-coef.aux.reftime);
    end
    v = [];
end

%%--------------------------------------------------------- 
function ut_diagnfigs(coef,indPE,t,u,v,usnrc,vsnrc,e)
% UT_DIAGNFIGS()
% create two diagnostic figures
% see report for a description of their contents
% UTide v1p0 9/2011 d.codiga@gso.uri.edu

if coef.aux.opt.twodim
    % first figure
    figure;
    % description
    subplot(4,1,1);
    text(0.2,0.6,coef.aux.rundescr,'fontsize',8);
    set(gca,'visible','off');
    % u time sequences
    subplot(4,1,2);
    plot(t,u,'r.:');
    hold;
    plot(t,usnrc,'b.:');
    plot(t,real(e),'g.:');
    set(gca,'xlim',[min(t) max(t)]);
    datetick('x','keepticks','keeplimits');
    ylabel('[raw input units]');
    title(sprintf(['UTide Diagnostic Plot #1:   u    Red= raw input;'...
        ' Blue= reconstructed fit (SNR>=%.1f); Green= residual'],...
        coef.aux.opt.diagnminsnr));
    grid on;
    % v time sequences
    subplot(4,1,3);
    plot(t,v,'r.:');
    hold;
    plot(t,vsnrc,'b.:');
    plot(t,imag(e),'g.:');
    set(gca,'xlim',[min(t) max(t)]);
    datetick('x','keepticks','keeplimits');
    ylabel('[raw input units]');
    title(sprintf(['v    Red= raw input; Blue= reconstructed fit'...
        ' (SNR>=%.1f); Green= residual'],coef.aux.opt.diagnminsnr));
    grid on;
    % all constituents, increasing frequency, rel SNR thresh
    subplot(4,1,4);
    [~,ind] = sort(coef.aux.frq,'ascend');
    SNRnum = coef.Lsmaj(ind).^2 + coef.Lsmin(ind).^2;
    SNRden = (coef.Lsmaj_ci(ind)/1.96).^2 + (coef.Lsmin_ci(ind)/1.96).^2;
    thresh = coef.aux.opt.diagnminsnr*SNRden;
    semilogy(thresh,'g');
    hold;
    ymin = 10.^(floor(log10(min([SNRnum; thresh])))-1);
    name = coef.name(ind);
    ind = find(SNRnum>=thresh);
    for i = 1:length(ind)
        line([ind(i) ind(i)],[ymin SNRnum(ind(i))],...
            'marker','*','color','r');
        text(ind(i),SNRnum(ind(i)),['  ' name(ind(i))],'rotation',60,...
            'color','r','verticalalign','bottom');
    end
    ind = find(SNRnum<thresh);
    for i = 1:length(ind)
        line([ind(i) ind(i)],[ymin SNRnum(ind(i))],...
            'marker','*','color','b');
        text(ind(i),SNRnum(ind(i)),['  ' name(ind(i))],'rotation',60,...
            'color','b','verticalalign','bottom');
    end
    set(gca,'ylim',[ymin max(get(gca,'ylim'))],'xlim',[0 length(SNRnum)+1],...
        'xticklabel',repmat(blanks(size(get(gca,'xticklabel'),2)),...
        size(get(gca,'xticklabel'),1),1));
    ylabel('L_{smaj}^{2} + L_{smin}^{2}');
    title(sprintf(['All constits, by increasing frequency; '...
        'Red= SNR>=%.1f, Blue= SNR<%.1f; Green= SNR thresh = %.1f('...
        '\\sigma_{L_{smaj}}^{2}+\\sigma_{L_{smin}}^{2})'],...
        coef.aux.opt.diagnminsnr,coef.aux.opt.diagnminsnr,...
        coef.aux.opt.diagnminsnr));
    % second figure
    figure;
    indPE = indPE(coef.diagn.SNR>=coef.aux.opt.diagnminsnr);
    % Lsmaj
    subplot(4,1,1);
    errorbar(coef.Lsmaj(indPE),coef.Lsmaj_ci(indPE),'.');
    grid on;
    set(gca,'xtick',1:length(indPE),'xticklabel',coef.name(indPE),...
        'xlim',[0 length(indPE)+1]);
    ylabel('[raw input units]');
    title(sprintf(['UTide Diagnostic Plot #2: Constits with SNR>%.1f, in '...
        'decreasing-PE order // Major axis L_{smaj} with 95%%'...
        ' confidence intervals'],coef.aux.opt.diagnminsnr));
    % Lsmin
    subplot(4,1,2);
    errorbar(coef.Lsmin(indPE),coef.Lsmin_ci(indPE),'.');
    grid on;
    set(gca,'xtick',1:length(indPE),'xticklabel',coef.name(indPE),...
        'xlim',[0 length(indPE)+1]);
    ylabel('[raw input units]');
    title('Minor axis L_{smin} with 95% confidence intervals');
    % orientation angle theta
    subplot(4,1,3);
    errorbar(coef.theta(indPE),coef.theta_ci(indPE),'.');
    grid on;
    set(gca,'ylim',[-45 225],'ytick',-45:45:225,'xlim',[0 length(indPE)+1],...
        'xtick',1:length(indPE),'xticklabel',coef.name(indPE));
    ylabel('[degrees]');
    title('Orientation angle theta with 95% confidence intervals');
    % phase lag g
    subplot(4,1,4);
    errorbar(coef.g(indPE),coef.g_ci(indPE),'.');
    grid on;
    set(gca,'ylim',[-90 450],'ytick',-90:90:450,'xlim',[0 length(indPE)+1],...
        'xtick',1:length(indPE),'xticklabel',coef.name(indPE));
    ylabel('[degrees]');
    title('Phase lag g with 95% confidence intervals');
    figure(1);
else % one-dimensional case
    % first figure
    figure;
    % description
    subplot(3,1,1);
    text(0.2,0.6,coef.aux.rundescr,'fontsize',8);
    set(gca,'visible','off');
    % time sequences
    subplot(3,1,2);
    plot(t,u,'r.:');
    hold;
    plot(t,usnrc,'b.:');
    plot(t,real(e),'g.:');
    set(gca,'xlim',[min(t) max(t)]);
    datetick('x','keepticks','keeplimits');
    ylabel('[raw input units]');
    title(sprintf(['UTide Diagnostic Plot #1:  Red= raw input;'...
        ' Blue= reconstructed fit (SNR>=%.1f); Green= residual'],...
        coef.aux.opt.diagnminsnr));
    grid on;
    % all constituents, increasing frequency, rel SNR thresh
    subplot(3,1,3);
    [~,ind] = sort(coef.aux.frq,'ascend');
    SNRnum = coef.A(ind).^2;
    SNRden = (coef.A_ci(ind)/1.96).^2;
    thresh = coef.aux.opt.diagnminsnr*SNRden;
    semilogy(thresh,'g');
    hold;
    ymin = 10.^(floor(log10(min([SNRnum; thresh])))-1);
    name = coef.name(ind);
    ind = find(SNRnum>=thresh);
    for i = 1:length(ind)
        line([ind(i) ind(i)],[ymin SNRnum(ind(i))],...
            'marker','*','color','r');
        text(ind(i),SNRnum(ind(i)),['  ' name(ind(i))],'rotation',60,...
            'color','r','verticalalign','bottom');
    end
    ind = find(SNRnum<thresh);
    for i = 1:length(ind)
        line([ind(i) ind(i)],[ymin SNRnum(ind(i))],...
            'marker','*','color','b');
        text(ind(i),SNRnum(ind(i)),['  ' name(ind(i))],'rotation',60,...
            'color','b','verticalalign','bottom');
    end
    set(gca,'ylim',[ymin max(get(gca,'ylim'))],'xlim',[0 length(SNRnum)+1],...
        'xticklabel',repmat(blanks(size(get(gca,'xticklabel'),2)),...
        size(get(gca,'xticklabel'),1),1));
    ylabel('A^{2}');
    title(sprintf(['All constits, by increasing frequency; '...
        'Red= SNR>=%.1f, Blue= SNR<%.1f; '...
        'Green= SNR thresh = %.1f \\sigma_{A}^{2} '],...
        coef.aux.opt.diagnminsnr,coef.aux.opt.diagnminsnr,...
        coef.aux.opt.diagnminsnr));
    % second figure
    figure;
    indPE = indPE(coef.diagn.SNR>=coef.aux.opt.diagnminsnr);
    % amplitude
    subplot(2,1,1);
    errorbar(coef.A(indPE),coef.A_ci(indPE),'.');
    grid on;
    set(gca,'xtick',1:length(indPE),'xticklabel',coef.name(indPE),...
        'xlim',[0 length(indPE)+1]);
    ylabel('[raw input units]');
    title(sprintf(['UTide Diagnostic Plot #2: Constits w/ SNR>%.1f, decreasing-PE order '...
        ' // Amplitude with 95%% confidence intervals'],...
        coef.aux.opt.diagnminsnr));
    % phase lag g
    subplot(2,1,2);
    errorbar(coef.g(indPE),coef.g_ci(indPE),'.');
    grid on;
    set(gca,'ylim',[-90 450],'ytick',-90:90:450,'xlim',[0 length(indPE)+1],...
        'xtick',1:length(indPE),'xticklabel',coef.name(indPE));
    ylabel('[degrees]');
    title('Phase lag g with 95% confidence intervals');
    figure(1);
end

%%--------------------------------------------------------- 
function coef = ut_finish(coef,nNR,nR,nI,elor,cnstit)
% UT_FINISH()
% 1 nan-fill coefficients & conf-ints if IRLS did not converge
%       (only portion that uses inputs other than coef)
% 2 create coef.results char array
% 3 order fields of coef
% 4 display coef.results, coef.aux.rundescr, and/or coef.diagn.table
%    as specified by opt.RunTimeDisp
% UTide v1p0 9/2011 d.codiga@gso.uri.edu

%% IRLS did not converge, so nan-fill coefficients and conf-ints
if ~isfield(coef,'g')
    if coef.aux.opt.twodim
        coef.Lsmaj = nan*ones(size(coef.name));
        coef.Lsmaj_ci = coef.Lsmaj;
        coef.Lsmin = coef.Lsmaj;
        coef.Lsmin_ci = coef.Lsmaj;
        coef.theta = coef.Lsmaj;
        coef.theta_ci = coef.Lsmaj;
        coef.g = coef.Lsmaj;
        coef.g_ci = coef.Lsmaj;
        coef.umean = NaN;
        coef.vmean = NaN;
        if ~coef.aux.opt.notrend
            coef.uslope = NaN;
            coef.vslope = NaN;
        end
    else
        coef.A = nan*ones(size(coef.name));
        coef.A_ci = coef.A;
        coef.g = coef.A;
        coef.g_ci = coef.A;
        coef.mean = NaN;
        if ~coef.aux.opt.notrend
            coef.slope = NaN;
        end
    end
    if ~coef.aux.opt.nodiagn
        nallc = nNR+nR+nI;
        nc = nNR+nR;
        coef.diagn.name = coef.name;
        coef.diagn.PE = nan*ones(nallc,1);
        coef.diagn.SNR = coef.diagn.PE;
        coef.diagn.TVraw = nan;
        coef.diagn.TVallc = nan;
        coef.diagn.TVsnrc = nan;
        coef.diagn.PTVallc = nan;
        coef.diagn.PTVsnrc = nan;
        coef.diagn.SNRallc = nan;
        coef.diagn.K = nan;
        coef.diagn.lo.name = cellstr(repmat('none',nallc,1));
        coef.diagn.hi.name = coef.diagn.lo.name;
        coef.diagn.lo.RR = nan*ones(nallc,1);
        coef.diagn.hi.RR = coef.diagn.lo.RR;
        [~,ind] = sort(coef.aux.frq(1:nc),'ascend');
        for i = 1:nc-1
            c1 = ind(i);
            c2 = ind(i+1);
            coef.diagn.lo.name(i+1) = coef.name(c1);
            coef.diagn.hi.name(i) = coef.name(c2);
            coef.diagn.lo.RR(i+1) = 24*elor*(coef.aux.frq(c2) - ...
                coef.aux.frq(c1))/coef.aux.opt.rmin;
        end
        coef.diagn.hi.RR(1:nc-1) = coef.diagn.lo.RR(2:nc);
        coef.diagn.lo.RNM = nan*ones(nallc,1);
        coef.diagn.hi.RNM = nan*ones(nallc,1);
        coef.diagn.lo.CorMx = nan*ones(nallc,1);
        coef.diagn.hi.CorMx = nan*ones(nallc,1);
    end
    minsnr = coef.aux.opt.diagnminsnr;
    tbl = cell(nallc+6+sign(nI)*(nI+1),1);
    tbl{1} = 'UTide Summary Diagnostics Table:';
    tbl{2} = sprintf(['   Rmin= %#-9.3g   MinSNR= %#-9.3g   '...
        '(* SNR >= MinSNR)'],coef.aux.opt.rmin,minsnr);
    tbl{3} = sprintf('   K= %#-9.3g   SNRallc= %#-9.3g',...
        coef.diagn.K,coef.diagn.SNRallc);
    tbl{4} = sprintf(['   TVallc= %#-9.3g    TVsnrc= %#-9.3g   '...
        'TVraw= %#-9.3g'],coef.diagn.TVallc,coef.diagn.TVsnrc,...
        coef.diagn.TVraw);
    tbl{5} = sprintf('   PTVallc= %5.1f%%   PTVsnrc= %5.1f%%',...
        coef.diagn.PTVallc,coef.diagn.PTVsnrc);
    tbl{6} = [' NAME    PE       SNR loNAME   loRR    loRNM  '...
        'loCorMx hiNAME   hiRR    hiRNM  hiCorMx'];
    main = [blanks(nallc)' char(coef.diagn.name)];
    main = [main reshape(sprintf('%6.2f%%',coef.diagn.PE),7,nallc)'];
    main = [main reshape(sprintf(' %#8.2g',coef.diagn.SNR),9,nallc)'];
    if nI
        coef.diagn.lo.name(end-nI+1:end) = ...
            cellstr(repmat('(I) ',nI,1));
        coef.diagn.hi.name(end-nI+1:end) = ...
            cellstr(repmat('(I) ',nI,1));
    end
    lo = char(coef.diagn.lo.name);
    lo = [lo reshape(sprintf(' %#8.2g',coef.diagn.lo.RR),9,nallc)'];
    lo = [lo reshape(sprintf(' %#8.2g',coef.diagn.lo.RNM),9,nallc)'];
    lo = [lo reshape(sprintf(' %#8.2g',coef.diagn.lo.CorMx),9,nallc)'];
    hi = char(coef.diagn.hi.name);
    hi = [hi reshape(sprintf(' %#8.2g',coef.diagn.hi.RR),9,nallc)'];
    hi = [hi reshape(sprintf(' %#8.2g',coef.diagn.hi.RNM),9,nallc)'];
    hi = [hi reshape(sprintf(' %#8.2g',coef.diagn.hi.CorMx),9,nallc)'];
    tbl(7:7+nallc-1) = cellstr([main repmat(' ',nallc,1) lo ...
        repmat(' ',nallc,1) hi]);
    if nI
        tbl{end-nI} = '  Inferred constituents (I):';
        cntr =0;
        for i=1:nR
            for j = 1:cnstit.R{i}.nI
                cntr = cntr+1;
                if coef.aux.opt.twodim
                    tbl{end-nI+cntr} = sprintf(['    %s inferred from'...
                        ' reference %s (r+,r- = %#.3g,%#.3g; ',...
                        'zta+,zta- = %#.3g,%#.3g deg)'],...
                        deblank(cnstit.R{i}.I.name{j}),...
                        deblank(cnstit.R{i}.name),...
                        abs(cnstit.R{i}.I.Rp(j)),...
                        abs(cnstit.R{i}.I.Rm(j)),...
                        (180/pi)*angle(cnstit.R{i}.I.Rp(j)),...
                        -(180/pi)*angle(cnstit.R{i}.I.Rm(j)));
                else
                    tbl{end-nI+cntr} = sprintf(['    %s inferred from'...
                        ' reference %s (r = %#.3g zta = %#.3g deg)'],...
                        deblank(cnstit.R{i}.I.name{j}),...
                        deblank(cnstit.R{i}.name),...
                        abs(cnstit.R{i}.I.Rp(j)),...
                        (180/pi)*angle(cnstit.R{i}.I.Rp(j)));
                end
            end
        end
    end
    coef.diagn.table = char(reshape(strrep(tbl(:)','NaN','---')',...
        size(tbl,1),size(tbl,2)));
end
%% create coef.results
coef.results = sprintf('UTide Results:');
if ~isnan(coef.g(1))
    nallc = length(coef.g);
    if coef.aux.opt.twodim
        coef.results = char(coef.results,[' Cnstit  Lsmaj  '...
            'Lsmaj_ci     Lsmin  Lsmin_ci     Theta  Theta_ci    '...
            '     g      g_ci']);
        main = [coef.Lsmaj coef.Lsmaj_ci coef.Lsmin coef.Lsmin_ci ...
            coef.theta coef.theta_ci coef.g coef.g_ci];
        for i = 1:nallc
            coef.results = char(coef.results,sprintf(['%4s %#9.3g '...
                '%#9.3g %#9.3g %#9.3g' ...
                ' %#9.3g %#9.3g %#9.3g %#9.3g'],coef.name{i},main(i,:)));
        end
        coef.results = char(coef.results,sprintf([' Mean u,v '...
            '(u/v units) = %#9.3g, %#9.3g'],coef.umean,coef.vmean));
        if ~coef.aux.opt.notrend
            coef.results = char(coef.results,sprintf([' Trend u,v slope '...
                '(u/v units per day) = %#9.3g, %#9.3g'],...
                coef.uslope,coef.vslope));            
        end
    else
        coef.results = char(coef.results,[' Cnstit      A'...
            '      A_ci         g      g_ci']);
        main = [coef.A coef.A_ci coef.g coef.g_ci];
        for i = 1:nallc
            coef.results = char(coef.results,sprintf(['%4s %#9.3g '...
                '%#9.3g %#9.3g %#9.3g'],coef.name{i},main(i,:)));
        end
        coef.results = char(coef.results,sprintf([' Mean (input '...
                'units) = %#9.3g' ],coef.mean));
        if ~coef.aux.opt.notrend
            coef.results = char(coef.results,sprintf([' Trend slope '...
                '(input units per day) = %#9.3g '],coef.slope));
        end
    end
else
    coef.results = char(coef.results,[' IRLS solution reached'...
        ' iteration limit without converging; '],...
        '    coefficients and confidence intervals set to NaN.');
end
%% order fields
if coef.aux.opt.twodim
    fldord = {'name'; 'Lsmaj'; 'Lsmaj_ci';...
        'Lsmin'; 'Lsmin_ci'; 'theta'; 'theta_ci'; 'g'; 'g_ci';...
        'umean'; 'vmean';};
else
    fldord = {'name'; 'A'; 'A_ci'; 'g'; 'g_ci';...
        'mean';};
end
if ~coef.aux.opt.notrend
    if coef.aux.opt.twodim
        fldord{end+1} = 'uslope';
        fldord{end+1} = 'vslope';
    else
        fldord{end+1} = 'slope';
    end
end
fldord{end+1} = 'results';
fldord{end+1} = 'aux';
if ~coef.aux.opt.nodiagn
    fldord{end+1} = 'diagn';
end
coef = orderfields(coef,fldord);
coef.aux = orderfields(coef.aux,{'rundescr';'opt';'frq';'lind';'lat';...
    'reftime';});
%% runtime display
if ~isequal(coef.aux.opt.runtimedisp,'nnn')
    fprintf('\n');
end
if isequal(coef.aux.opt.runtimedisp(1),'y')
    disp(char(coef.results,' '));
end
if isequal(coef.aux.opt.runtimedisp(2),'y')
    disp(char(coef.aux.rundescr,' '));
end
if isequal(coef.aux.opt.runtimedisp(3),'y')
    if ~coef.aux.opt.nodiagn && ~isnan(coef.g(1))
        disp(char(coef.diagn.table,' '));
    end
end

