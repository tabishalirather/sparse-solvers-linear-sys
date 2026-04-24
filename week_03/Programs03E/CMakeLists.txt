
# Require a recent version of CMake.
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

# Set project name.
project("Itsolver" C) 

# Enable or disable verbose output when building targets.
# Default behaviour is disabled.
set(CMAKE_VERBOSE_MAKEFILE FALSE)

# Turn on compiler warnings.
set(CMAKE_C_FLAGS "-Wall -Wextra -Wno-unused-function")

# Define USE_BLAS if cmake was invoked with -DUSE_BLAS=1 on the command line.
if(USE_BLAS)
  add_definitions("-DUSE_BLAS")
endif()

# Link to cairo if specified.
if(USE_CAIRO)
  link_libraries("cairo")
endif()

# Include directory where the compiler will search for header files.
include_directories(Library)

# Add subdirectory containing library sources.
add_subdirectory(Library)

# Enable testing meaning that CMake will actually generate tests when the
# add_test command is used. Calls to add_test would have no effect otherwise.
# The entire test suite can be executed using
#   make test
# on the command line.
# Tests are considered passed when the program returns with exit code zero.
# Each test program is run in the directory where the executable file is
# located in the build directory.
enable_testing()

# Copy matrix market files to build directory to ensure test program with
# hardcoded paths to these files will find them when run in the build
# directory.
#file(COPY "gre_512.mtx" DESTINATION ".")
#file(COPY "e05r0500.mtx" DESTINATION ".")
#file(COPY "e05r0500_rhs1.mtx" DESTINATION ".")
#file(COPY "orsirr_1.mtx" DESTINATION ".")

# Function to be invoked in CMakeLists.txt files in subdirectories.
# The definition of this function has to precede the calls to add_subdirectory
# to assure that the function will be visible in subdirectories.
# Adds executable targets coming from source files named
#   testprog_*.c
# residing within the directory from which the function is invoked.
function(itsolver_add_testprogs)
  # Get a list of all C source files.
  file(GLOB sources *.c)

  # Collect C source files corresponding to test programs.
  # Initialize list to be empty since it might already contain
  # files that were inserted in a higher-level CMakeList.txt.
  unset(testprog_sources)
  foreach(source ${sources})
    if(${source} MATCHES testprog_.*\\.c$)
      set(testprog_sources ${testprog_sources} ${source})
    endif()
  endforeach()

  # Generate executable test program targets.
  foreach(testprog_source ${testprog_sources}) 
    # Remove .c suffix from testprog source file. 
    string(REGEX MATCH testprog_.*\\.c$ testprog ${testprog_source})
    string(REGEX REPLACE \\.c$ "" testprog ${testprog})
    add_executable(${testprog} ${testprog_source})

    # Link with itsolver and C math library.
    target_link_libraries(${testprog} itsolver m)

    # Link with BLAS only when specified.
    if(USE_BLAS)
      target_link_libraries(${testprog} openblas)
    endif()

    # Add program to the list of tests to be run automatically.
    add_test(NAME ${testprog} COMMAND ${testprog})
  endforeach()
endfunction(itsolver_add_testprogs)

# Call the previously defined function to add test programs from the current directory.
itsolver_add_testprogs()

# Add subdirectories, notable after the preceding function definition.
# Get a list of subdirectories.
file(GLOB subdirs Examples*)
# Note that we have to loop over all subdirectories.
# Giving the list of subdirectories to add_subdirectory directly
# does not work because additional arguments to this function have
# a special significance.
foreach(subdir ${subdirs})
  add_subdirectory(${subdir})
endforeach()
