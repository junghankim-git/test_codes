
import os


var = [0.0, 0.0]


def test_func(a):
#    a = 0
#    b = [1, 2]
    global var
    var = a
    return var



def main():
    aa = test_func([1.0,1.0])
    #print('aa = {}, {}'.format(aa[0],aa[1]))
    print 'aa = ', aa, var
    bb = test_func([2.0,2.0])
    print 'bb = ', bb, var
    
    #print('aa = {}, bb = {}'.format(aa,bb))
    print 'aa = ', aa, ', bb = ', bb, ', var = ', var
    
    
    
    '''
    
    aa = -1
    bb = [0, 0]
    
    test_func(aa, bb)
    
    print 'aa = ', aa
    print 'bb = ', bb
    '''

    return 0


err = main()


