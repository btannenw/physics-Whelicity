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
                #print systVar.split(',')[0], systVar.split(',')[2], systVar.split(',')[4].split('\n')[0]
                name = systVar.split(',')[0]
                if name != 'PartonShower' and name != 'MCgenerator' and name != 'ColorReconnection' and name != 'UnderlyingEvent':
                    os.system("""root -l -q -b 'calcStatError.C("{0}","{1}","{2}", "{3}")'""".format( systVar.split(',')[0], systVar.split(',')[2], systVar.split(',')[4].split('\n')[0], "" ))
                else:
                    passBool = 'true'
                    os.system("""root -l -q -b 'calcStatError.C("{0}","{1}","{2}", "{3}")'""".format( systVar.split(',')[0], systVar.split(',')[2], systVar.split(',')[4].split('\n')[0], passBool ))
            else:
                print "$$$$$$$$$$$$$$ NON 5"
                print systVar.split(',')[0], systVar.split(',')[2], systVar.split(',')[4].split('\n')[0]
                os.system("""root -l -q -b 'calcStatError.C("{0}","{1}","{2}", "{3}")'""".format( systVar.split(',')[0], systVar.split(',')[0], systVar.split(',')[0], "" ))
