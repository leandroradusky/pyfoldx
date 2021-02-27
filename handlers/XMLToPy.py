'''
XML2Py - XML to Python de-serialization

This code transforms an XML document into a Python data structure

Usage:
    deserializer = XML2Py()
    python_object = deserializer.parse( xml_string )
    print xml_string
    print python_object
'''

from lxml import etree

class XML2Py():

    def __init__( self ):

        self._parser = parser = etree.XMLParser( remove_blank_text=True ) #@UnusedVariable
        self._root = None  # root of etree structure
        self.data = None   # where we store the processed Python structure

    def parse( self, xmlString ):
        '''
        processes XML string into Python data structure
        '''
        self._root = etree.fromstring( xmlString.encode("ascii") , self._parser )
        self.data = self._parseXMLRoot()
        return self.data

    def tostring( self ):
        '''
        creates a string representation using our etree object
        '''
        if self._root != None:
            return etree.tostring( self._root )

    def _parseXMLRoot( self ):
        '''
        starts processing, takes care of first level idisyncrasies
        '''
        childDict = self._parseXMLNode( self._root )
        return { self._root.tag : childDict["children"] }

    def _parseXMLNode( self, element ):
        '''
        rest of the processing
        '''
        childContainer = None # either Dict or List

        # process any tag attributes
        # if we have attributes then the child container is a Dict
        #   otherwise a List
        if element.items():
            childContainer = {}
            childContainer.update( dict( element.items() ) )
        else:
            childContainer = {}


        if isinstance( childContainer, list ) and element.text:
            # tag with no attributes and one that contains text
            childContainer.append( element.text )
        
        else:
            # tag might have children, let's process them
            for child_elem in element.getchildren():

                childDict = self._parseXMLNode( child_elem )

                # let's store our child based on container type
                #
                if isinstance( childContainer, dict ):
                    # these children are lone tag entities ( eg, 'copyright' )
                    if childDict["tag"] in childContainer:
                        if isinstance( childContainer[childDict["tag"]] ,list ):
                            childContainer.update( { childDict["tag"] : childContainer[childDict["tag"]]+[childDict["children"]] } )
                        else:
                            childContainer.update( { childDict["tag"] : [childContainer[childDict["tag"]]]+[childDict["children"]] } )
                    else:
                        childContainer.update( { childDict["tag"] : childDict["children"] } )

                else:
                    # these children are repeated tag entities ( eg, 'format' )
                    childContainer.append( childDict["children"] )
            if element.text != None:
                if element.text.strip() != "":
                    childContainer['text'] = element.text.replace('\n','')

        return { "tag":element.tag, "children": childContainer }