import sys
import subprocess
import glob


def convertArc(arc):
    startVertex = arc.split(",")[0].split("(")[1]
    endVertex = arc.split(",")[1].split(")")[0]
    return startVertex, endVertex

def preprocess(graph):
    preprocessedGraph = []
    for i in range(0, len(graph)):
        if graph[i][0] == "(": preprocessedGraph.append(graph[i].split(' ')[0][:-1])
    return preprocessedGraph

def printGraph(inputName, outputName):
    inputFile = open(inputName, "r")
    outputFile = open(outputName, "w")
    graph = inputFile.readlines()
    graph = preprocess(graph)
    outputFile.write("digraph G {\n")
    for arc in graph:
        startVertex, endVertex = convertArc(arc)
        outputFile.write("  " + startVertex + " -> " + endVertex + ";\n")
    outputFile.write("}")


def outputPath(inputFilePath):
    return "./Images/" + inputFilePath.split("/")[2][:-4] + ".gv"

if __name__ == '__main__':
    fileList = glob.glob("./Pre-Test/*")
    for inputFile in fileList:
        outputFile = outputPath(inputFile)
        printGraph(inputFile, outputFile)
    outputFileList = glob.glob("./Images/*")
    for outputFile in outputFileList:
        subprocess.call("dot -Tpng " + outputFile + " -o " + outputFile[:-3] + ".png", shell = True)
