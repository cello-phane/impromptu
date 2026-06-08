#!/bin/bash
gcc -Wall -O2 -o impromptu.exe impromptu.c -I"C:/Program Files (x86)/SDL2/include/SDL2" -L"C:/Program Files (x86)/SDL2/lib" -lsdl2 -lmingw32 -lSDL2main -lSDL2 -mwindows

#build test_radicaltrig.c
#gcc -Wall -std=c99 -O2 test_radicaltrig.c radicaltrig.c -lm -o test_radicaltrig -I"C:/Program Files (x86)/SDL2/include/SDL2" -L"C:/Program Files (x86)/SDL2/lib" -lSDL2 -lm
