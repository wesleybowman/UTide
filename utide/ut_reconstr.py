from ut_reconstr1 import ut_reconstr1


def ut_reconstr(tin, coef, **opts):

    u, v = ut_reconstr1(tin, coef, **opts)

    return u, v
