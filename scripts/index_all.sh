# indexes all FASTA files in a given dir

# takes a directory and samtools path

# check if the user provided a directory
if [ -z "$1" ]; then
    echo "Please provide a directory"
    exit 1
fi

# check if the user provided a samtools path
if [ -z "$2" ]; then
    echo "Please provide a samtools path"
    exit 1
fi

# check if the directory exists
if [ ! -d "$1" ]; then
    echo "Directory does not exist"
    exit 1
fi

# get the directory
dir=$1

for file in $dir/*.fasta; do
    if [ ! -f "$file.fai" ]; then
      echo "Indexing $file"
        ${2} faidx $file
    else
      echo "Index exists for file $file"
    fi

    # dict
    if [ ! -f "${file}.dict" ]; then
      echo "Creating dict for $file"
        ${2} dict $file -o ${file}.dict -u .
    else
      echo "Dict exists for file $file"
    fi
done
