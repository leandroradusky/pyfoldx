'''
Created on Mar 25, 2013

@author: leandro
@summary: this class allows user to interact with files
'''

import gzip
import os.path, os
from src.handlers.ErrorHandler import ErrorHandler
from src.handlers.SystemHandler import SystemHandler


class FileHandler(object):
    '''
    classdocs
    '''
    
    def __init__(self):
        '''
        Constructor
        '''
    
    @staticmethod
    def fileExists(path):
        try:
            with open(path): return True
        except IOError:
            return False
    
    @staticmethod
    # @param path: the path of the file to be readed
    # @return: a list of string with all the lines, [] if error 
    def getLines(path):
        f = None
        try:
            f = open(path, 'r')
            lines = f.readlines()
            f.close()
            ret = []
            for l in lines:
                ret += [l.strip()]
            return ret
        except Exception as inst:
            if f != None:
                f.close()
            ErrorHandler.printError("handlefile.getLines", {"path": path}, inst)
            return []
    
    @staticmethod
    # @param path: the path of the file to be readed
    # @return: a list of string with all the lines, [] if error 
    def getGZLines(path):
        f = None
        try:
            f = gzip.open(path, 'r')
            lines = f.readlines()
            f.close()
            ret = []
            for l in lines:
                ret += [l.strip()+"\n"]
            return ret
        except Exception as inst:
            if f != None:
                f.close()
            ErrorHandler.printError("handlefile.getLines", {"path": path}, inst)
            return []
    
    @staticmethod
    # @param path: the path of the file to be readed
    # @return: a string with the file 
    def getString(path):
        try:
            return "\n".join(FileHandler.getLines(path))
        except Exception as inst:
            ErrorHandler.printError("handlefile.getString", {"path": path,}, inst)
            return ""

    @staticmethod
    # @param path: the path of the file where the line will be appended
    # @param line: a string to append in a new line 
    # @return: True is success, False if don't
    def appendLine(path, line):
        f = None
        try:
            FileHandler.ensureDir(path)
            f = open(path, 'a')
            f.write(line+"\n")
            f.close()
            return True
        except Exception as inst:
            if f != None:
                f.close()
            ErrorHandler.printError("handlefile.appendLine", {"path": path, "line":line}, inst)
            return False
    
    @staticmethod
    # @param path: the path of the file where the line will be appended
    # @param lines: a list of strings to append in a new line 
    # @return: True is success, False if don't
    def appendLines(path, lines):
        f = None
        try:
            FileHandler.ensureDir(path)
            f = open(path, 'a')
            for line in lines:
                f.write(line+"\n")
            f.close()
            return True
        except Exception as inst:
            if f != None:
                f.close()
            ErrorHandler.printError("handlefile.appendLines", {"path": path}, inst)
            return False
    
    @staticmethod
    # @param path: the path of the file where the line will be written
    # @param line: a string to write in a new line 
    # @return: True is success, False if don't
    def writeLine(path, line, addReturn = "\n"):
        f = None
        try:
            FileHandler.ensureDir(path)
            f = open(path, 'w')
            f.write(line+addReturn)
            f.close()
            return True
        except Exception as inst:
            ErrorHandler.printError("handlefile.writeLine", {"path": path, "line":line}, inst)
            if f != None:
                f.close()
            return False
    
    @staticmethod
    # @param path: the path of the file where the line will be written
    # @param line: a string to write in a new line 
    # @return: True is success, False if don't
    def writeLines(path, lines):
        f = None
        try:
            FileHandler.ensureDir(path)
            f = open(path, 'w')
            for line in lines:
                if len(line)== 0 or line[-1] == "\n":
                    f.write(line)
                else:
                    f.write(line+"\n")
            f.close()
            os.chmod(path,0o755)
            return True
        except Exception as inst:
            if f != None:
                f.close()
            ErrorHandler.printError("handlefile.writeLines", {"path": path, "line":lines}, inst)
            return False
    
    @staticmethod
    # @param path: the dir to make if not exist
    # @summary: make dir if not exists 
    def ensureDir(path):
        d = os.path.dirname(path)
        if not os.path.exists(d):
            os.makedirs(d)
            SystemHandler.executeCommand("chmod -R 777 %s" % d)
        return d
    
    @staticmethod
    # @param path: the dir to make if not exist
    # @summary: make dir if not exists 
    def exists(path):
        return os.path.exists(path)
    
    