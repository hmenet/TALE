
import numpy as np

#some auxiliary functions used by both sampling and reconciliation likelihood computation


#functions for logspace addition and substraction
def log_add(a,b):
    m=min(a,b)
    M=max(a,b)
    return M + np.log(1+np.exp(m-M))

#a>b
def log_minus(a,b):
    if a >= b:
        return a+np.log(1-np.exp(b-a))

def log_add_list(l):
    M=max(l)
    s=0
    for a in l:
        s+=(np.exp(a-M))
    return M+np.log(s)


#to check if leaves have multiple matches
def is_mult_match(c_match):
    return type(next(iter(c_match.values()))) == list
#next return the next element of an iterator
#iter enable dict, that is iterable, to have a next function, so O(1).

#test if intersection of lists l1 and l2 is not empty
def inter_list(l1,l2):
    s1=set(l1)
    return(any(x in s1 for x in l2))
#any return true if any element of an iterable is true