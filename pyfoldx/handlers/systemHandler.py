'''
Created on Mar 25, 2013

@author: leandro
@summary: this class allows user to interact with system command line
'''

import os

class SystemHandler(object):
    '''
    Statatic methods to work with the system in a friendly manner.
    '''

    def __init__(self):
        '''
        Constructor: this is a class of static methods.
        '''
        
    @staticmethod
    def executeCommand(command):
        '''
        Execute a command.
        
        :param command: the command to be executed.
        '''
        try:
            os.system(command)
        except Exception as inst:
            pass
        
    @staticmethod
    # @param command: the command to be executed
    # @return: a list of the lines returned by the command, [] if error
    def getCommandResult(command, strip_lines = False):
        '''
        Execute a command and return its output.
        
        :param command: the command to be executed.
        :param strip_lines: split the output in a list of strings, one for each line.
        :return: a list of the lines returned by the command, [] if error
        '''
        try:
            import warnings
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                lines = os.popen(command).readlines()
                if strip_lines:
                    lines = map(lambda x: x.strip(), lines)
        except Exception as inst:
            pass
            return []
        return lines
