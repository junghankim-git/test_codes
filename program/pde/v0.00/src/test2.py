
import os



def test_func(a, b):
    a = 0
    b[0] = 1
    b[1] = 2



def main():
    
    aa = -1
    bb = [0, 0]
    
    test_func(aa, bb)
    
    print 'aa = ', aa
    print 'bb = ', bb

    return 0


err = main()


