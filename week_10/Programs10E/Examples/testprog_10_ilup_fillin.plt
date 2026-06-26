set terminal pdf

set output 'testprog_10_ilup_fillin.pdf'

set xlabel 'p'
set ylabel 'nnz(A) / n^2'

set key outside box

set title 'Comparison of ILU(p) for different orderings'

plot 'testprog_10_ilup_fillin.txt' using 1:2 title 'no reordering' with linespoints pointtype 5 pointsize .25,\
'testprog_10_ilup_fillin.txt' using 1:3 title 'Cuthill McKee' with linespoints pointtype 5 pointsize .25,\
'testprog_10_ilup_fillin.txt' using 1:4 title 'Nested Dissection' with linespoints pointtype 5 pointsize .25,\
'testprog_10_ilup_fillin.txt' using 1:5 title 'Independent Set' with linespoints pointtype 5 pointsize .25
