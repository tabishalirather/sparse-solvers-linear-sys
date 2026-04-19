
# Initialise list of library source files.
file(GLOB library_sources *.c)

# Add library subdirectories where additional library source files will be
# added to the list.

# Itsolver library target.
add_library(itsolver ${library_sources})

