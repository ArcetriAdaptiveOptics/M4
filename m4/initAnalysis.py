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
from m4.analyzers import timehistory as th


###_____________________________________________________
line = ["",
        "Using the IPython console for M4 data analysis",
        "                                 _______      ",
        "M       M         44            /\     /\     ",
        "M M   M M        4 4           /  \   /  \    ",
        "M   M   M       4  4          /    \ /    \   ",
        "M       M      4 4 4 4       ||-----X-----||  ",
        "M       M          4          \    / \    /   ",
        "M       M          4           \  /   \  /    ",
        "                                \/_____\/     "
   ]
for l in line:
    print(l)
