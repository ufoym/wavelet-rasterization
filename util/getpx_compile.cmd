@echo off
g++ -c -o getpx.obj getpx.cpp
g++ -shared -o getpx.dll getpx.obj
del getpx.obj
@pause