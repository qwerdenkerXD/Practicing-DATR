name: Franz-Eric Sill
student number: 2139315

how to run the script?
    1. You need python3 installed on your device so in terminal you get a version number 3.x.y of python if you type in: python3 -V
       On Ubuntu it's preinstalled, on other distributions you may have to do it manually by your package manager or sth. others. There are plenty of guides in the internet of how to do it.
       During development Python 3.8.10 was used.
    2. In terminal run: python3 -m pip install matplotlib numpy  # these are the necessary modules for the script
       During development matplotlib version 3.5.2 and numpy version 1.23.0 were used.
    3. Open the terminal in the folder where the "Sill_Franz_Eric_2139315.zip" was extracted (the folder of this readme.txt, where the SL.py lies), so the command "pwd" should return the path of this folder.
    4. Make sure that in the folder are no files named "heatmap.png" and "trees.png", because they will be overwritten, and run the following command: python3 SL.py

output: "trees.png" with the tree plots only and "heatmap.png" with the heatmap and attached trees