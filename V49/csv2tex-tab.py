#Usage:
#csv2textable
#<pathToFile>:  Input the path to your input file
#<caption>:     The caption you'd like to have for the table
#<label>:       The label you'd like to have for the table. This is also the name for the output file
#<format>:      The table format. For expample cccc for a for four columns with centered entries
#<delimiter>:   The delimiter of the csv-file. If you don't give this argument, the standard value is set to ";"
#
#The csv-file should be formated like this:
#headline1;headline2;headline3
#entry11;entry12;entry13
#entry21;entry22;entry33
#etc
#
#You should have an equal amount of colums in every row. The rows should not end with your delimiter.
#
#You will need to have a directory named "content" relative to the location of this python-file.
#This script probably won't work on Windows systems. Install a linux distribution.
#
#csv_read
#<pathToFile>:  Input the path to your input file
#<delimiter>:   The delimiter of the csv-file. If you don't give this argument, the standard value is set to ";"

def csv_read(pathToFile, delimiter=";"):
    with open(pathToFile, "r") as f:
        content = []
        for line in f:
            content.append((line.rstrip()).split(delimiter))
    return content

def csv2textable(pathToFile, caption, label, format, delimiter=";"):
    content = csv_read(pathToFile, delimiter)
    
    tableHead = "\\begin{table}[!htp]\n" + "\\centering\n" + "\\caption{" + caption + "}\n" + "\\label{tab:" + label + "}\n"
    tabularHead = "\\begin{tabular}{" + format + "}\n" + "\\toprule\n"

    fileEnd = "\\end{tabular}\n" + "\\end{table}"

    with open("content/" + label + ".tex", "w") as f:
        f.write(tableHead)
        f.write(tabularHead)
        for i in range(len(content)):
            for j in range(len(content[i])):
                if(i == 0):
                    f.write("{" + content[i][j] + "} ")
                    if(j == (len(content[i]) - 1)):
                        f.write("\\\\\n")
                        f.write("\\midrule\n")
                    else:
                        f.write("& ")
                else:
                    f.write(content[i][j])
                    if(j == (len(content[i]) - 1)):
                        f.write(" \\\\\n")
                    else:
                        f.write(" & ")
        
        f.write("\\bottomrule\n")
        f.write(fileEnd)


file="data/diff.csv"
caption="Gemessene Spannungen in Abh??ngigkeit von $\tau$ f??r die Bestimmung des Diffusionskoeffizienten."
label="diff"
form="c c"
csv2textable(file, caption, label, form)