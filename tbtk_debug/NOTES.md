# Compilation

Example compile command

```
/usr/bin/c++  \
-DTBTK_MATPLOTLIB_DO_NOT_FORCE_AGG -DTBTK_USE_OPEN_MP -DTBTK_VERSION_GIT_HASH=\"\" \
-DTBTK_VERSION_MAJOR=2 -DTBTK_VERSION_MINOR=6 \
-DTBTK_VERSION_PATCH=4 -DTBTK_WRAP_PRIMITIVE_TYPES=1 \
-I/usr/include/python3.8 \
-I/home/efefer/WORKS/TBTK_test/MyApp/include \
-I/home/efefer/mysoftwares/tbtk-2.6.4/include \
-fopenmp -std=c++11 -Wall -O3 -o main.o -c main.cpp
```

# Linking

```
/usr/bin/c++   -fopenmp -std=c++11 -Wall -O3  -rdynamic main.o  -o build/main.x \
-L/home/efefer/mysoftwares/tbtk-2.6.4/lib/TBTK/CMake/..  \
-Wl,-rpath,/home/efefer/mysoftwares/tbtk-2.6.4/lib/TBTK/CMake/.. \
-lTBTK -lfftw3 -lgmp -lgmpxx -lpython3.8 -llapack -lblas
```