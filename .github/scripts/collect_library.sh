# library from .rmd
grep 'library(' FigureYa*/*.Rmd > library.txt
sed 's/ //g' library.txt | sed 's/"//g' | sed 's/\// /g' | sed 's/(/ /g' | sed 's/)/ /g' | awk '{print $1,$3}' | sort -k 2| sort -k 1 | uniq > FigureYa_library.txt
cut -d' ' -f2 FigureYa_library.txt | sort | uniq > library2download.txt

# github from install_dependencies.R 
grep 'https://github.com/' FigureYa*/install_dependencies.R | sed -E 's/.*(https:\/\/.*)/\1/' > github_library_R.txt
# github from .Rmd
grep 'https://github.com/' FigureYa*/*.Rmd | grep -oE 'https?://[^ ]+' | sed 's/[)>].*//g' | sed 's/[一-龥]//g' | sort | uniq > github_library_rmd.txt
cat github_library_R.txt github_library_rmd.txt | cut -d'/' -f1-5 > github_library.txt

# other source from install_dependencies.R 
grep 'https://' FigureYa*/install_dependencies.R > other_library_R.txt

grep -vE 'https://bioconductor.org/|https://cloud.r-project.org/|https://github.com/|https://cran.r-project.org/|https://cloud.r-project.org' other_library_R.txt | grep -oE 'https?://[^ ]+' | sed 's/[)>].*//g' | sed 's/[一-龥]//g' > filtered_file.txt
