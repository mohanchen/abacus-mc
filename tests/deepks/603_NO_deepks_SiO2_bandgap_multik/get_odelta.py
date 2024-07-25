import numpy
a=numpy.load('OUT.autotest/otot.npy')
b=numpy.load('OUT.autotest/obase.npy')
print((a-b)[0][0])
