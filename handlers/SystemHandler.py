'''
Created on Mar 25, 2013

@author: leandro
@summary: this class allows user to interact with system command line
'''

from src.handlers.ErrorHandler import ErrorHandler
import os
import argparse

class SystemHandler(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''
        
    @staticmethod
    # @param command: the command to be executed
    # @return: a list of the lines returned by the command, [] if error
    def executeCommand(command):
        try:
            os.system(command)
        except Exception as inst:
            ErrorHandler.printError(__name__, {"command": command}, inst)
        
    @staticmethod
    # @param command: the command to be executed
    # @return: a list of the lines returned by the command, [] if error
    def getCommandResult(command, stripLines = False):
        try:
            lines = os.popen(command).readlines()
            if stripLines:
                lines = map(lambda x: x.strip(), lines)
        except Exception as inst:
            ErrorHandler.printError(__name__, {"command": command}, inst)
            return []
        return lines
    
    @staticmethod
    # @param params: list of tuples with (name, desc) of each command line argument
    # @return: object with each param as a member
    def getCommandLineParameters(params):
        try:
            parser = argparse.ArgumentParser()
            for param in params:
                parser.add_argument(param[0], type=str,
                                    help=param[1])
            
            args = parser.parse_args()
            return args
        except Exception as inst:
            ErrorHandler.printError(__name__, {"params": params}, inst)
            return None
        