#!/usr/bin/python

from math import ceil, floor
from operator import itemgetter
import csv
import sys


current_word = None
current_count = 0
word = None

# input comes from STDIN
for line in sys.stdin:
    # remove leading and trailing whitespace
    line = line.strip()

    # parse the input we got from mapper.py
    all = line.split('\t', 1)
    count = 1
    word = all[0]

    try:
        count = int(count)
    except ValueError:
        # count was not a number, so silently
        # ignore/discard this line
        continue

    # this IF-switch only works because Hadoop sorts map output
    # by key (here: word) before it is passed to the reducer
    if current_word == word:
        current_count += count
    else:
        if current_word:
            now = current_word
            want = now.split(',')
            x_a = float(want[0]) - 0.1
            y_a = float(want[1]) - 0.1
            # write result to STDOUT
            print '%s\t%s\t%s\t%s\t%s' % (x_a, want[0] , y_a, want[1], current_count)
        current_count = count
        current_word = word

# do not forget to output the last word if needed!
#if current_word == word:
#	print '%s\t%s' % (current_word, current_count)
