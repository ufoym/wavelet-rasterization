@echo off
g++ -c -o solver.obj solver.cpp
g++ -shared -o solver.dll solver.obj
del solver.obj