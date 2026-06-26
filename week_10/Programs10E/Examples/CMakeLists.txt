# Add testprograms
itsolver_add_testprogs()

# Copy gnuplot scripts to the build directory
file(GLOB gnuplot_scripts *.plt)
foreach(script ${gnuplot_scripts})
    file(COPY ${script} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endforeach(script ${gnuplot_scripts})