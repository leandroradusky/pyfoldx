'''
Created on Jul 12, 2013
@author: leandro
'''
from pyfoldx.handlers.fileHandler import FileHandler
from pyfoldx.handlers.systemHandler import SystemHandler
from datetime import datetime
PROXY = ""

class URLRetrieveHandler(object):
    
    @staticmethod
    def RetrieveFileLines(url):
        a = datetime.utcnow()
        jobid=str(a.year)+str(a.month)+str(a.day)+str(a.hour)+str(a.minute)+str(a.second)+str(a.microsecond)
        URLRetrieveHandler.RetrieveFileToDisk(url, "/tmp/", jobid)
        lines = FileHandler.getLines("/tmp/"+jobid)
        SystemHandler.executeCommand("rm /tmp/"+jobid)
        return lines
        
    @staticmethod
    def RetrieveFileToDisk(url,path, filename=""):
        dir = FileHandler.ensureDir(path)
        command = 'cd %s; wget "%s" %s 2> /dev/null' % \
        (dir,url, "" if filename == "" else " -O  %s " % filename)
        SystemHandler.getCommandResult(command)



