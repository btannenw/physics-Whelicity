#!/usr/bin/python

import os, subprocess, sys

systMin = int(sys.argv[1])
systMax = int(sys.argv[2])

systList = open ('systConfigFull.txt','r')

i=0
for systVar in systList:
    if '#' not in systVar:
        i=i+1
	if i >= systMin and i <= systMax:
            #systVar = systVar.split(',')[0]
            print len(systVar.split(',')), systVar.split(',')
            if len(systVar.split(',')) == 5:
                print "$$$$$$$$$$$$$$ YEAH 5"
                #print systVar.split(',')[0], systVar.split(',')[2], systVar.split(',')[4].split('\n')[0]
                os.system("""root -l -q -b 'calcStatError.C("{0}","{1}","{2}")'""".format( systVar.split(',')[0], systVar.split(',')[2], systVar.split(',')[4].split('\n')[0] ))
            else:
                print "$$$$$$$$$$$$$$ NON 5"
                print systVar.split(',')[0], systVar.split(',')[2], systVar.split(',')[4].split('\n')[0]
                os.system("""root -l -q -b 'calcStatError.C("{0}","{1}","{2}")'""".format( systVar.split(',')[0], systVar.split(',')[0], systVar.split(',')[0] ))
