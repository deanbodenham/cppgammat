testexecfile="test_exec"

g++ -std=c++11 -isystem ~/build/googletest/googletest/include/ -pthread \
    ./main.cpp  \
    ../gammatest/src/utils.cpp \
    ./test1.cpp  \
    ~/build/googletest/googletest/libgtest.a -o $testexecfile

./$testexecfile
