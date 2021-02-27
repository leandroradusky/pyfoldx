'''
Created on Nov 2, 2020
@author: lradusky
'''

from src.structure.Structure import Structure

if __name__ == '__main__':
    
    code = "1exp"
    #st = Structure(code, "/data/cif/%s.cif" % code)
    #st = Structure(code, "/home/lradusky/Downloads/%s.pdb" % code)
    st = Structure(code)
    print( st.toPdb() )
    
    print( "done" )