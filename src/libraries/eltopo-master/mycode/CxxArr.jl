using Cxx

r = [1,2,3]
x = convert(cxxt"std::vector< int >", r)


cxxinclude("vector")
cxx" std::vector<int> a = std::vector<int> (5,6); "
c = icxx" &a[0]; "
cSize = icxx" a.size(); "
unsafe_wrap(Array,c,cSize,true)

### Mixed example


r = Int32[1,2,3]
cxxinclude("vector")
rc = convert(cxxt"std::vector< int >", r)
rp = icxx" &$rc[0]; "
rSize = icxx" $rc.size(); "
rr = unsafe_wrap(Array,rp,rSize,true)

### So in principle I can avoid of writting the fucking loops.
### I will get a proper exceptions
### I do not need to hardcode the size of array.
### I do not need to compile and write the wrapper outside the julia

### The next step so is a to learn using an external code




using Cxx

# Importing shared library and header file
const path_to_lib = pwd()
addHeaderDir(path_to_lib, kind=C_System)
Libdl.dlopen(path_to_lib * "/libArrayMaker.so", Libdl.RTLD_GLOBAL)
cxxinclude("ArrayMaker.h")

# Creating class object
maker = @cxxnew ArrayMaker(5, 2.0)

arr = @cxx maker->fillArr()

## So this indeed works
unsafe_wrap(Array,arr,5)
#pointer_to_array(arr, 5)

### I could use Cxx now and figure out how to build a shared library of ElTopo
### Then with that I could proceed on learning to use BinaryBuilder


#### How teh wrapper could look like
const path_to_lib = pwd()
addHeaderDir(path_to_lib, kind=C_System)
Libdl.dlopen(path_to_lib * "/eltopo.so", Libdl.RTLD_GLOBAL)
#cxxinclude
