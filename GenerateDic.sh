#!/bin/bash
#### For Cern-root
### This Generates a Dictianry file with .so related to LinkDef.h, then you need to use gSystem->Load("libPDvector.so") to load it to your code
####
rootcling -v4 -f PDvector.cxx  -rmf libPDvector.rootmap -rml libPDvector.so  LinkDef.h 
g++ -fPIC -shared -o libPDvector.so PDvector.cxx `root-config --cflags --libs --glibs `


