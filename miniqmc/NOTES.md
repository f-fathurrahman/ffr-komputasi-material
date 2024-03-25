Example compilation command:
```
g++ -I ../../src -I ../../catch \
-fopenmp -fomit-frame-pointer -finline-limit=1000 \
-fstrict-aliasing -funroll-all-loops -D__forceinline=inline \
-Wno-deprecated -march=native -O3 -DNDEBUG -ffast-math   -std=c++11  \
-DADD_ -DH5_USE_16_API -DHAVE_CONFIG_H -Drestrict=__restrict__ -c \
ParticleSet.cpp
```

