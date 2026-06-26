
# Initialise list of library source files.
file(GLOB library_sources *.c)

# Add library subdirectories where additional library source files will be
# added to the list.
add_subdirectory(CSR)
add_subdirectory(LinearSolver)
add_subdirectory(Models)
add_subdirectory(Permutations)
add_subdirectory(Plot)
add_subdirectory(KrylovSolver)
add_subdirectory(Factorizations)

# Itsolver library target.
add_library(itsolver ${library_sources})

