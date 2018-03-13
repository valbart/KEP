# Script for auto testing performances of the post-ex evaluation of the pre-testing model
# Executable, configuration file, solution and instance files HAVE to be in the same directory
# This script HAS to be launch from the same directory (relative path are used)

import sys
import glob
import subprocess

# Name of the executable to be tested

# No ordering nor graph splitting during post-ex evaluation
RKEPOriginal = "RobustKEPSolverOriginal"
# Ordering of the vertices (k.e.a order at start then order by in/out degree)
RKEPOrderingOnly = "RobustKEPSolverCurrentBest"
# Ordering the vertices after splitting the graph using articulation points
RKEPSpliting = "RobustKEPSolverTest"

# This function modify the config file in order to launch the post-ex evaluation
# on the instance instanceFileName with solution file solutionFileName
# The original configuration file must have its solver option set to 5 (exact expected
# transplant) as we do not check for this option here.
# Name of the config file HAS to be "0Config.txt"
def modifyConfigFile(instanceFileName, solutionFileName):
    configFile = open("0Config.txt", "r")
    configLines = configFile.readlines()
    for i in range(0, len(configLines)):
        if configLines[i].split("=")[0] == "Input File ":
            configLines[i] = "Input File = " + instanceFileName + "\n"
        if configLines[i].split("=")[0] == "Solution Input File ":
            configLines[i] = "Solution Input File = " + solutionFileName + "\n"
    configFile.close()
    configFile = open("0Config.txt", "w")
    configFile.truncate()
    for line in configLines: configFile.write(line)
    configFile.close()

# Run the executable exec on the 0Config.txt configuration file.
# Input and solution file must have been set before caling this function
# Output the computation time in "TmpOutput.txt" and return the "real time"
# computation in number of seconds as a string
def runKEP(execName):
    r = subprocess.call( "(time ./" + execName + " 0Config.txt) 2> TmpOutput.txt", shell = True)
    if r == 1: return "Bug"
    tmpFile = open("TmpOutput.txt", "r")
    lines = tmpFile.readlines()
    seconds = int(lines[1].split("\t")[1].split("\n")[0].split("m")[1].split("s")[0].split(",")[0])
    minutes = int(lines[1].split("\t")[1].split("\n")[0].split("m")[0])
    tmpFile.close()
    return (str(60*minutes + seconds))

# Return the instance file that corresponds to the given solution file. To do so,
# we identify the instance number x (Said-x) and the failure probability that figure
# in the name of the solution file.
def indentifyInstance(solutionFileName):
    instanceNb = solutionFileName[-5]
    if (instanceNb == '0'): instanceNb = "10" # Instance number goes from 1 to 10
    pr = solutionFileName.split("-")[3]
    return "V-50-" + pr + "-Said-" + instanceNb + ".txt"


if __name__ == '__main__':
    solutionFileList = glob.glob("V-PT-*.txt") # Only list solution files
    outputFile = open("OutputAutoTest.txt", "w")
    outputFile.write("# Output format: Solution file name, time original, time ordering, time spliting\n\n")
    for solutionFileName in solutionFileList:
        instanceFileName = indentifyInstance(solutionFileName)
        modifyConfigFile(instanceFileName, solutionFileName)
        timeOriginal = runKEP(RKEPOriginal)
        timeOrdering = runKEP(RKEPOrderingOnly)
        timeSpliting = runKEP(RKEPSpliting)
        outputFile.write(solutionFileName + "," + timeOriginal + "," + timeOrdering + "," + timeSpliting + "\n")
    outputFile.close()
