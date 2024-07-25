import numpy
a=numpy.load('OUT.autotest/htot.npy')
b=numpy.load('OUT.autotest/hbase.npy')
print(numpy.sum(a-b))
