def partition_generator(n,k,l):
    if k < 1:
        raise StopIteration
    if k == 1:
        if n >= l:
            yield (n,)
        raise StopIteration
    for i in range(l,n+1):
        for result in partition_generator(n-i,k-1,i):
            yield (i,)+result
