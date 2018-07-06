basnam=`basename $1 .tex`
lualatex -shell-escape $1
#bibtex $basnam
#lualatex -shell-escape $1
#lualatex -shell-escape $1
