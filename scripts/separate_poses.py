#!/usr/bin/env python

# C. Langini - modified 13/06/2017

#import math
#import re
import sys
import string
import getopt
from itertools import islice

class Pose(object):
    def __init__(self, name, poseNum, cluster, fr_nu, totEn):
        self.name = str(name)
        self.poseNum = int(poseNum)
        self.cluster = int(cluster)
        self.fr_nu = int(fr_nu)
        self.totEn = float(totEn)
    def __repr__(self):
        return "Name: %30s Pose: %8d Clu: %8d Fr_nu: %8d TotEn: %8.2f" %(
                        self.name, self.poseNum, self.cluster, self.fr_nu, self.totEn)
    def __eq__(self,other):
        if isinstance(other, Pose): # check if dealing with istances of same object
            return (self.name == other.name) and (self.poseNum == other.poseNum) and \
                    (self.cluster == other.cluster) and (self.totEn == other.totEn) and \
                    (self.fr_nu == other.fr_nu)
        return NotImplemented

def dump_mol2(fIn,header, pose, poseRank):
    #print "Found pose: " + str(pose)
    #outBaseName = "_pose_"
    outFn = pose.name + "_" + str(poseRank) + ".mol2"
    #print outFn
    with open(outFn,'w') as fOut:
        fOut.write("@<TRIPOS>MOLECULE\n")
        fOut.write(pose.name + "\n")
        for headerLine in header:
            fOut.write(headerLine)
        for line in fIn:
            if line.strip() == "@<TRIPOS>MOLECULE":
                inNextMol = True
                return inNextMol
            fOut.write(line)
            #fOut.write(line)
            #fOut.write("\n")

def main(argv):
    inFn = ''
    tabFn = ''
    npos = ''
    try:
          opts, args = getopt.getopt(argv,"hi:t:n:",["imol2=","table=","nposes="])
    except getopt.GetoptError:
        print 'separate_poses.py -i <pproc mol2> -t <seed summary table> -n <number of poses>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'separate_poses.py reads the first n lines from a SEED output table (please remove header line!) and\n',\
                    'extracts the corresponding poses from the concatenated mol2 saving them to ',\
                    'separate mol2 files.\nUsage:'
            print '\tpython separate_poses.py -i <pproc mol2> -t <seed summary table> -n <number of poses>'
            print '\tpython separate_poses.py --imol2=<pproc mol2> --table=<seed summary table> --nposes=<number of poses>'
            sys.exit()
        elif opt in ("-i", "--imol2"):
            inFn = arg
            print "Input mol2: %s" % inFn
        elif opt in ("-t", "--table"):
            tabFn = arg
            print "SEED summary table: %s" % tabFn
        elif opt in ("-n", "--nposes"):
            npos = int(arg)
            print "Number of poses to extract: %d" % npos

    if (inFn != '') and (tabFn != '') and (npos != ''):
        #print "Extracting poses"
        poseList = []
        poseRank = range(1,npos+1)
        with open(tabFn, 'r') as f:
            for line in islice(f,npos):
                pose_split = line.strip().split()
                curPose = Pose(pose_split[0], pose_split[1], pose_split[2], pose_split[3], pose_split[4])
                poseList.append(curPose)
    # OLDER VERSION:
        # with open(tabFn, 'r') as f:
        #     seedTab = list(islice(f, npos))
        # for line in seedTab:
        #     #poseList = [Pose() for count in xrange(4)]
        #     pose_split = line.strip().split()
        #     curPose = Pose(pose_split[0],pose_split[1],pose_split[2],pose_split[3])
        #     poseList.append(curPose)
        # #print poseList

        found_all = False
        nfound = 0
        #molStart = None # pointer to molecule start
        inNextMol = False
        with open(inFn,'r') as fIn:
            while (not found_all):
                if (not inNextMol):
                    for line in fIn:
                        if line.strip() == "@<TRIPOS>MOLECULE":
                            break
                    else:
                        break
                inNextMol = False
                for line in islice(fIn,1):
                    pass
                curName = line.strip()
                #print curName
                header = list(islice(fIn, 5))
                #print header
                #line = next(iter(fIn))
                line = header[-1]
                #print line
                pose_split = line.strip().split()
                curPose = Pose(curName,pose_split[1],pose_split[3],pose_split[5], pose_split[7])
                #print "current: " + str(curPose)
                for pose_idx, pose in enumerate(poseList):
                    #print pose
                    if (pose == curPose):
                        nfound = nfound + 1
                        print "Found pose: " + str(pose)
                        inNextMol = dump_mol2(fIn,header,pose,poseRank[pose_idx])
                        del poseList[pose_idx]
                        del poseRank[pose_idx]
                        break
            # OLDER VERSION:
                # for pose in poseList:
                #     #print pose
                #     if (pose == curPose):
                #         nfound = nfound + 1
                #         print "Found pose: " + str(pose)
                #         inNextMol = dump_mol2(fIn,header,pose,nfound)
                #         poseList.remove(pose)
                #         break
                if (nfound == npos):
                    found_all = True
        if found_all:
            print "All poses extracted"
        else:
            print "EOF reached after extracting " + str(nfound) + " poses."
            print "Missed poses: "
            for pose in poseList:
                print pose



    else:
        print "Missing parameters."
        print 'separate_poses.py reads the first n lines from a SEED output table (please remove header line!) and\n',\
                'extracts the corresponding poses from the concatenated mol2 saving them to ',\
                'separate mol2 files.\nUsage:'
        print '\tpython separate_poses.py -i <pproc mol2> -t <seed summary table> -n <number of poses>'
        print '\tpython separate_poses.py --imol2=<pproc mol2> --table=<seed summary table> --nposes=<number of poses>'



if __name__ == '__main__':
    main(sys.argv[1:])
