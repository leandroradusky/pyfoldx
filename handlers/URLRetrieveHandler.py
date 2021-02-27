'''
Created on Jul 12, 2013
@author: leandro
'''
from src.handlers.FileHandler import FileHandler
from src.handlers.SystemHandler import SystemHandler

import urllib

PROXY = ""

class URLRetrieveHandler(object):
    
    @staticmethod
    def RetrieveFileLines(url):
        SystemHandler.executeCommand("rm /tmp/download")
        URLRetrieveHandler.RetrieveFileToDisk(url, "/tmp/", "download")
        return FileHandler.getLines("/tmp/download")
        
    @staticmethod
    def RetrieveFileToDisk(url,path, filename=""):
        dir = FileHandler.ensureDir(path)
        command = "cd %s; /usr/local/bin/curl --ciphers DEFAULT@SECLEVEL=1 %s -o %s 2> /dev/null" % (dir,url,filename)
        SystemHandler.executeCommand(command)
        



