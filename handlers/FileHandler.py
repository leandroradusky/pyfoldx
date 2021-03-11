'''
Created on Mar 25, 2013

@author: leandro
@summary: this class allows user to interact with files
'''

import gzip, os
from pyfoldx.handlers.SystemHandler import SystemHandler


class FileHandler(object):
    '''
    Statatic methods to work with files in a friendly manner.
    '''
    
    def __init__(self):
        '''
        Constructor: this is a class of static methods.
        '''
    
    @staticmethod
    def fileExists(path):
        '''
        Check if a file exists.
        
        :param path: The file or folder to check if exists.
        :return: True if the file or folder exists. 
        '''
        try:
            with open(path): return True
        except IOError:
            return False
    
    @staticmethod
    def getLines(path):
        '''
        Get lines of a file.
        
        :param path: the path of the file to be read.
        :return: a list of string with all the lines, [] if error 
        '''
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
            return []
    
    @staticmethod
    def getGZLines(path):
        '''
        Get lines of a gzipped file.
        
        :param path: the path of the file to be read.
        :return: a list of string with all the lines, [] if error 
        '''
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
            return []
    
    @staticmethod
    def getString(path):
        '''
        Get string with the contents of a file.
        
        :param path: the path of the file to be read.
        :return: a string with all the contents of a file, "" if error 
        '''
        try:
            return "\n".join(FileHandler.getLines(path))
        except Exception as inst:
            return ""

    @staticmethod
    def appendLine(path, line):
        '''
        Append a line to a file.
        
        :param path: the path of the file where the line will be appended
        :param line: a string to append in a new line 
        :return: True is success, False if don't
        '''
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
            return False
    
    @staticmethod
    # @param path: the path of the file where the line will be appended
    # @param lines: a list of strings to append in a new line 
    # @return: True is success, False if don't
    def appendLines(path, lines):
        '''
        Append a list of lines to a file.
        
        :param path: the path of the file where the line will be appended
        :param lines: a list of strings to append in a new line 
        :return: True is success, False if don't
        '''
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
            return False
    
    @staticmethod
    def writeLine(path, line, addReturn = "\n"):
        '''
        Write a line into a file.
        
        :param path: the path of the file where the line will be appended
        :param line: a string to append in a new line 
        :return: True is success, False if don't
        '''
        f = None
        try:
            FileHandler.ensureDir(path)
            f = open(path, 'w')
            f.write(line+addReturn)
            f.close()
            return True
        except Exception as inst:
            if f != None:
                f.close()
            return False
    
    @staticmethod
    def writeLines(path, lines):
        '''
        Write a list of lines to a file.
        
        :param path: the path of the file where the line will be appended
        :param lines: a list of strings to append in a new line 
        :return: True is success, False if don't
        '''
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
            return False
    
    @staticmethod
    def ensureDir(path):
        '''
        Ensure the existance of a directory.
        
        :param path: the dir to make if not exist
        :summary: make dir if not exists 
        '''
        d = os.path.dirname(path)
        if not os.path.exists(d):
            os.makedirs(d)
            SystemHandler.executeCommand("chmod -R 777 %s" % d)
        return d
    
    @staticmethod
    def exists(path):
        '''
        Check the existance of a file or directory.
        
        :param path: the dir or file to check if exist
        :return: True if exists, false otherwise
        '''
        return os.path.exists(path)
    
    
    
    


