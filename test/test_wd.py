from __future__ import unicode_literals

import doctest
import unittest
import wd


class TestFunc(unittest.TestCase):
    def test_get_Nh(self):
        self.assertEqual(wd.get_Nh(0), 0) 
        print('test_get_Nh2 OK')
    
    def test_get_Nh2(self):
        self.assertEqual(wd.get_Nh(100), 10.0) 
        print('test_get_Nh2 OK')
    
    def test_aprox(self):        
        self.assertEqual(wd.aprox([17769.0,16865.743295019158],[54552,64552.0],0,60000),17276.90574712644)  
        print('test_aprox OK')
       
        
if __name__ == '__main__':
    unittest.main()      
