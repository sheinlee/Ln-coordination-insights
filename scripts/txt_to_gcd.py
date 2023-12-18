import os

files = os.listdir('.')
for filename in files:
    portion = os.path.splitext(filename)
    if portion[1] == '.txt':
        newname = portion[0] +'.gcd'
        os.rename(filename,newname)