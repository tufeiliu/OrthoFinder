#! /bin/sh

read -p "Please provide the database dir: " database

if [ ! -d "$database" ]; then
    echo "Database not found: $database"
    exit 1
fi


read -p "Please provide the file path that contains the criteria: " filepath

if [ ! -f "$filepath" ]; then
    echo "File not found: $filepath"
    exit 1
fi
 
read -p "Please provide the distination data dir: " datadir


# Proteomes16.txt
while IFS= read -r line|| [ -n "$line" ]; do

    trimmed_line="${line%"${line##*[![:space:]]}"}"
    filename="${trimmed_line##*/}"

    for file in "$database"/*; do
        if [ -f "$file" ]; then

            database_filename=$(basename "$file")
            if [[ "$filename" == "$database_filename" ]]; then
                cp "$file" "$datadir"
                echo "copy $filename to $datadir"
            fi

        fi
    done
done < "$filepath"