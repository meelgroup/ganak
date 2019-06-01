import argparse
import os
import re
import sys
import time

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('DIMACSCNF', nargs='?', type=str, default="", help='input cnf file')
    parser.add_argument("-q", action='store_true', help="quiet mode", dest='q')
    parser.add_argument("-t ", type=int, default=5000, help="set time bound to T seconds", dest='t')
    parser.add_argument("-maxcachesize ", type=int, default=2000, help="set max cache size to CS MB", dest='cs')
    parser.add_argument("-noPP", action='store_true', help="turn off preprocessing", dest='noPP')
    parser.add_argument("-noCC", action='store_true', help="turn off component caching", dest='noCC')
    parser.add_argument("-noIBCP", action='store_true', help="turn off implicit BCP", dest='noIBCP')
    parser.add_argument("-noPCC", action='store_true', help="turn off probablistic component caching (PCC)", dest='PCC')
    parser.add_argument("-noCSVSADS", action='store_true', help="turn off CSVSADS branching heuristic", dest='CSVSADS')
    parser.add_argument("-noPC", action='store_true', help="turn off polarity cache heuristic (PC)", dest='PC')
    parser.add_argument("-noIS", help="turn off independent support heuristic (IS)", dest='IS', action='store_true')
    parser.add_argument("-noLSO", help="turn off learn and start over heuristic (LSO)", dest='NOLSO', action='store_true')
    parser.add_argument("-EDR", help="use exponentially decaying randomness (EDR)", dest='EDR', action='store_true')
    parser.add_argument("-lso ", type=int, default=5000, help="use learn and start over heuristic after LSO decisions", dest='LSO')
    parser.add_argument("-delta ", type=float, default=0.05, help="the confidence parameter", dest='DELTA')
    parser.add_argument("-seed", type=int, default=1000, help="seed for randomness", dest='seed')
    args = parser.parse_args()
    total_user_time = 0
    mis_calculated = False
    if not (args.DIMACSCNF):
        parser.error("Please provide the CNF formula file")
    if (args.CSVSADS and args.EDR):
        parser.error("CSVSADS and EDR cannot be used together ")
    cmd = "/usr/bin/time --verbose -o timeoutfile ./ganak -cs " + str(args.cs) + " -t " + str(args.t) + " "
    cmd += " -seed "+ str(args.seed) + " " 
    if (args.noCC):
        cmd += "-noCC "
    if (args.noIBCP):
        cmd += "-noIBCP "
    if (args.noPP):
        cmd += "-noPP "
    if (args.q):
        cmd += "-q "
    if (not args.PCC):
        cmd += "-pcc "
    if (not args.CSVSADS):
        cmd += "-act "
    if (not args.PC):
        cmd += "-pol random "
    else:
        cmd += "-pol default "
    if (not args.NOLSO):
        cmd += "-res " + str(args.LSO) + " "
    if (args.EDR):
        cmd += "-rand "  
    
    #setting the value of delta
    cmd += "-delta " + str(args.DELTA) + " "
    new_cmd = cmd
    if (not args.IS):
        cmd += "-maxdec " + "1000000" + " " + "500" + " "   #number of decisions and number of conflicts
    
    #setting the range of clhash family (= m*64)
    m = 1
    cmd += "-m " + str(m) + " "

    cmd += args.DIMACSCNF + " > outputfile 2>&1 "
    os.system(cmd)
    f = open("timeoutfile")
    text = f.read()
    f.close()
    total_user_time += float(re.findall(r"User time.*",text)[0].split(":")[1])

    f = open("outputfile")
    output  = f.read().strip().split("\n")[-1]
    f.close()
    if (output.startswith("Terminating")):
        mis_cmd = "/usr/bin/time --verbose  python2.7 MIS.py -timeout=100 -useInd=0 -output=mis.out " + args.DIMACSCNF + " > mis.timeout 2>&1"
        os.system(mis_cmd)
        f = open("mis.timeout")
        text = f.read()
        f.close()
        total_user_time += float(re.findall(r"User time.*",text)[0].split(":")[1])
        cmd = new_cmd
        cmd += "-m " + str(m) + " "
        cmd += "-is mis.out " 
        cmd += args.DIMACSCNF + " > outputfile 2>&1 "
        os.system(cmd)
        f = open("timeoutfile")
        text = f.read()
        f.close()
        total_user_time += float(re.findall(r"User time.*",text)[0].split(":")[1])
        f = open("outputfile")
        output = f.read().strip().split("\n")[-1]
        f.close()
        mis_calculated = True
    while (output == "-1"):
        m = 2 * m               #Double the range
        print("Change the hash range to 64x",m)
        cmd = new_cmd
        cmd += "-m " + str(m) + " "
        if(mis_calculated):
            cmd += "-is mis.out "
        elif (not args.IS):
            cmd += "-maxdec " + "1000000" + " " + "500" + " " #number of decisions and number of conflicts
        cmd += args.DIMACSCNF + " > outputfile 2>&1 "
        os.system(cmd)
        f = open("timeoutfile")
        text = f.read()
        f.close()
        total_user_time += float(re.findall(r"User time.*",text)[0].split(":")[1])
        f = open("outputfile")
        output = f.read().strip().split("\n")[-1]
        f.close()
        if (output.startswith("Terminating")):
            mis_cmd = "/usr/bin/time --verbose  python2.7 MIS.py -timeout=100 -useInd=0 -output=mis.out " + args.DIMACSCNF + " > mis.timeout 2>&1"
            os.system(mis_cmd)
            f = open("mis.timeout")
            text = f.read()
            f.close()
            total_user_time += float(re.findall(r"User time.*",text)[0].split(":")[1])
            cmd = new_cmd
            cmd += "-m " + str(m) + " "
            cmd += "-is mis.out " 
            cmd += args.DIMACSCNF + " > outputfile 2>&1 "
            os.system(cmd)
            f = open("timeoutfile")
            text = f.read()
            f.close()
            total_user_time += float(re.findall(r"User time.*",text)[0].split(":")[1])
            f = open("outputfile")
            output = f.read().strip().split("\n")[-1]
            f.close()
    #clear the extra files
    os.system("rm mis.out mis.timeout timeoutfile")
    os.system("cat outputfile")
    print("The total user time taken by GANAK is: ", total_user_time)

if __name__== "__main__":
    main()
