"""
Author(s)
---------
        - Pietro Ferraiuolo: written in 2024
        
Description
-----------
Script that initialize the ipython shell to be data analysis ready.
"""
from m4.utils import osutils
from m4.configuration import update_folder_paths as ufp
from m4.mini_OTT import timehistory as th




###_____________________________________________________
line = ["",
        " M4 data analysis ",
        "     _______      ",
        "    /\     /\     ",
        "   /  \   /  \    ",
        "  /    \ /    \   ",
        " ||-----X-----||  ",
        "  \    / \    /   ",
        "   \  /   \  /    ",
        "    \/_____\/     "
   ]

for l in line:
    print(l)
    