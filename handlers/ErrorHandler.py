'''
Created on Mar 25, 2013

@author: leandro
@summary: easy error handling
'''

import inspect

PRINT_ERRORS = False

class ErrorHandler(object):
    '''
    classdocs
    '''
    
    def __init__(self):
        '''
        Constructor
        '''
    
    @staticmethod
    # @param function:the function calling 
    # @param data: the function parameter values
    # @param exception: the error itself 
    def printError(function, data, exception):
        if PRINT_ERRORS:
            print( "-----------------------------------------" )
            print( "Error executing ", function )
            print( "with: " ) 
            for d in data:
                print( d, " = ", data[d] )
            print( "-----------------------------------------" )
            print( "Error Details:" )
            print( type(exception) )    # the exception instance
            print( exception.args )     # arguments stored in .args
            print( exception )          # __str__ allows args to printed directly
            print( "-----------------------------------------" )
            print( "-----------------------------------------" )
            print( "Error Traceback:" )
            for elem in inspect.stack(): 
                print( elem )
            print( "-----------------------------------------")
        