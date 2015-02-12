import numpy as np


def ut_rcninit(tin, **opts):

    t = tin[:]

    t[np.isnan(t)] = []
    # t(isnan(t)) = []
    opt = {}

    opt['cnstit'] = False
    opt['minsnr'] = 2
    opt['minpe'] = 0

    for key, item in opts.items():
        try:
            opt[key] = item
        except KeyError:
            print('ut_reconstr: unrecognized input: {0}'.format(key))

    # args = list(args)
    # args = [string.lower() for string in args]

    # Need an example of the args

#    while ~isempty(args)
#        switch(lower(args{1}))
#            case 'cnstit'
#                opt.cnstit = args{2};
#                args(1:2) = [];
#            case 'minsnr'
#                opt.minsnr = args{2};
#                args(1:2) = [];
#            case 'minpe'
#                opt.minpe = args{2};
#                args(1:2) = [];
#            otherwise
#                error(['ut_reconstr: unrecognized input: ' args{1}]);
#        end
#    end

    return t, opt
