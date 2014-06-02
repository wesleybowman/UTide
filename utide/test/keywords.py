def test(*args,**kwargs):
    print kwargs
    # print a
    print kwargs['a']
    opt={}
    opt['a']= 10
    opt['b']= 10
    opt['c']= 10
    print opt

    for i,v in kwargs.items():
        print i,v
        opt[i]=v
    test2(**kwargs)
    print opt

def test2(**kwargs):
    print 'hi'

def test3(a=10, b=10):
    test3.keywords=['a','b','c']
    print 't3'
    print a
    tt={}
    tt['a']=a
    tt['b']=b

    return tt



x=[]
y=[]
z=[]
test(x,y,z, a=2,b=3)

print test3.__defaults__
tt = test3(a=4, b=6)
print test3.__defaults__
print test3.__dict__
print test3.__repr__()
print 'HI'
test(**tt)
