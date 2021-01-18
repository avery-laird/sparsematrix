# Build 

### 1) Setup up python

Create a virtual environment (recommend)

Use pip to install the packages in `requirements.txt`

### 2) Create Makefile

```
mkdir build 
cd build
cmake --target bench
```

CMake will automatically run `gen_triang.py`, which generates the code required by `bench`.

Finally, make bench:

```
make bench
```

### 3) Adding Matrices

There are two options for including additional matrices into the build process.

The first option is to add the paths for the matrix and vector to `gen_triang.py`. This will generate appropriate code 
for the matrix into `optimized_triang_test<test-numer>.cpp`. This can then be included in `bench.cpp`, where the
function `int optsolve<test-number>(int, int*, int*, int*, double*&)` can be used.

The second option is to run the command:

```
python gen_triang.py <path/to/matrix/A> <path/to/vector/b> > <desired_file>.cpp 
```

And then use `<desired_file>.cpp` as required.

# Test

The file `tests.py` contains some useful tests:

```
pytest tests.py
```
