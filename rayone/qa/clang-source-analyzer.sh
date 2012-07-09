srcdir=../src
bindir=.
src=$srcdir/rayone.cpp
bin=$bindir/rayone
#cc=g++
scn=~/lab/checker-267/scan-build
cc=clang++
cmd="$scn $cc -o $bin $src\
	-std=c++0x -O4 -Wfatal-errors -Wall -Wextra -Werror -Wconversion -pedantic -pedantic-errors\
	-framework OpenGL\
	-framework GLUT\
	-I/System/Library/Frameworks/GLUT.framework/Versions/Current/Headers\
	-I/System/Library/Frameworks/OpenGL.framework/Versions/Current/Headers\
	-L/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries
"&&
echo $cmd&&
$cmd&>clang-source-analyzer.log&&
rm $bin&&
echo