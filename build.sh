#!/bin/bash
gcc -Wall -O2 -o impromptu.exe impromptu.c -I"C:/Program Files (x86)/SDL2/include/SDL2" -L"C:/Program Files (x86)/SDL2/lib" -lsdl2 -lmingw32 -lSDL2main -lSDL2 -mwindows
