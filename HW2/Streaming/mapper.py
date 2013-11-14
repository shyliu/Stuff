#!/usr/bin/python

from math import ceil, floor
import sys

# input comes from STDIN (standard input)

# input comes from STDIN (standard input)
for line in sys.stdin:
    line = line.strip()
    words = line.split()


    more = []
    for w in range(len(words)):
        new = float(words[w])
        new = ceil(new*10)/10
        more.append(new)
    C = map(str, more)
    character = C[0]+ ',' +C[1]
    print '%s\t%s' % (character, 1)




