set pm3d; set palette
set palette color
set pm3d map
set cbrange [:]
set xrange [-15:100]
set yrange [0:1]
unset ztics
set samples 101
set isosamples 20
set xtics 10
set ytics 0.1

set term jpeg 6
set output "gexp_3d.jpeg"

set palette model RGB

set palette defined 
set title "set palette defined"
splot 'f_SN.dat'
