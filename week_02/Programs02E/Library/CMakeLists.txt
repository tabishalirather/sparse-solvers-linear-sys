
# Initialise list of library source files.
file(GLOB library_sources *.c)

# Add library subdirectories where additional library source files will be
# added to the list.
add_subdirectory(CSR)
add_subdirectory(Models)
add_subdirectory(Permutations)

# Itsolver library target.
#add_library(itsolver ${library_sources}
#        Permutations/perm-not-efficient.c)
add_library(itsolver ${library_sources}
        CSR/csr_permute_code_rot.c)
